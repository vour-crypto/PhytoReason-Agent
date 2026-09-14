"""PhytoReason local research workbench (PySide6)."""
from __future__ import annotations

import html
import json
import os
import sys
from datetime import datetime, timedelta, timezone
from pathlib import Path

_PROJECT_ROOT = Path(__file__).resolve().parents[2]

try:
    from PySide6.QtCore import QObject, QThread, Qt, Signal
    from PySide6.QtGui import QFont
    from PySide6.QtWidgets import (QApplication, QCheckBox, QComboBox, QFileDialog,
        QFormLayout, QFrame, QGridLayout, QHBoxLayout, QLabel, QLineEdit, QListWidget,
        QListWidgetItem, QMainWindow, QMessageBox, QPushButton, QDockWidget,
        QScrollArea, QSplitter, QStackedWidget, QTextBrowser, QVBoxLayout, QWidget,
        QInputDialog)
except ImportError:  # pragma: no cover
    QApplication = None


def _parse_chips(evidence_json: str) -> list[str]:
    try:
        value = json.loads(evidence_json or "[]")
    except Exception:
        return []
    return [str(item) for item in value][:4] if isinstance(value, list) else []


_ENV_LLM_KEYS = {"api_key": "OPENAI_API_KEY", "base_url": "OPENAI_BASE_URL", "model": "LLM_MODEL"}


def _legacy_llm_settings(db) -> dict:
    """按 config.json > WorkbenchDB > env 的优先级解析 LLM 配置。

    Phase 6.4 Step 2：Web 与 Legacy 收口为单一事实源——config.json 优先；
    WorkbenchDB 仅作为历史设置的中间回退层；env 只兜底默认值。
    """
    from phyto_reason.config.llm_config import DEFAULTS, file_config
    file_cfg = file_config()
    resolved: dict[str, str] = {}
    for key in ("provider", "base_url", "model", "api_key"):
        if file_cfg.get(key):
            resolved[key] = str(file_cfg[key])
            continue
        db_value = db.get_setting(key, "")
        if db_value:
            resolved[key] = db_value
            continue
        env_name = _ENV_LLM_KEYS.get(key)
        env_value = os.getenv(env_name, "") if env_name else ""
        resolved[key] = env_value or str(DEFAULTS.get(key, ""))
    return resolved


if QApplication is not None:
    class ClickableCard(QFrame):
        clicked = Signal()
        def __init__(self, parent=None) -> None:
            super().__init__(parent); self.setObjectName("card"); self.setCursor(Qt.CursorShape.PointingHandCursor)
        def mousePressEvent(self, event) -> None:
            if event.button() == Qt.MouseButton.LeftButton: self.clicked.emit()
            super().mousePressEvent(event)

    class ConfidenceBar(QFrame):
        """6px confidence track with sage-to-tan fill."""
        def __init__(self, percent: float, parent=None) -> None:
            super().__init__(parent); self.percent = max(0.0, min(100.0, float(percent))); self.setFixedHeight(6); self.setStyleSheet("background:#e6e3dc;border-radius:3px;"); self._fill = QFrame(self); self._fill.setStyleSheet("background:qlineargradient(x1:0,y1:0,x2:1,y2:0,stop:0 #d4e0ce,stop:1 #e6cdb1);border-radius:3px;")
        def resizeEvent(self, event) -> None:
            self._fill.setGeometry(0, 0, max(3, int(self.width() * self.percent / 100)), 6); super().resizeEvent(event)

    class HypothesisCard(QFrame):
        RATING_CLASS = {"Plausible": "plausible", "Weak": "weak", "Insufficient": "insufficient"}
        card_clicked = Signal(str)
        def __init__(self, hypothesis_id: str, statement: str, rating: str, confidence: float | None, chips: list[str] | None = None, detailed: bool = False, parent=None) -> None:
            super().__init__(parent); self.hypothesis_id = hypothesis_id; self.setObjectName("card"); self.setCursor(Qt.CursorShape.PointingHandCursor); layout = QVBoxLayout(self); layout.setContentsMargins(18, 16, 18, 16); layout.setSpacing(7)
            pair = QLabel(f"🧬 {statement}"); pair.setObjectName("pairTitle"); pair.setWordWrap(True); layout.addWidget(pair)
            badge = QLabel(rating); badge.setObjectName("ratingBadge"); badge.setProperty("rating", self.RATING_CLASS.get(rating, "insufficient")); badge.setAlignment(Qt.AlignmentFlag.AlignCenter); badge.setFixedWidth(100); layout.addWidget(badge)
            value = float(confidence or 0); percent = value * 100 if value <= 1 else value; layout.addWidget(ConfidenceBar(percent)); note = QLabel(f"置信度 {percent:.0f}%"); note.setObjectName("muted"); layout.addWidget(note)
            if detailed:
                row = QHBoxLayout()
                for chip in chips or []: label = QLabel(chip); label.setObjectName("chip"); row.addWidget(label)
                row.addStretch(); layout.addLayout(row)
        def mousePressEvent(self, event) -> None:
            if event.button() == Qt.MouseButton.LeftButton: self.card_clicked.emit(self.hypothesis_id)
            super().mousePressEvent(event)

    class AgentWorker(QObject):
        event = Signal(dict); failed = Signal(str); finished = Signal()
        def __init__(self, agent, message: str, session_id: str) -> None: super().__init__(); self.agent, self.message, self.session_id = agent, message, session_id
        def execute(self) -> None:
            try:
                for event in self.agent.handle_stream(self.message, self.session_id): self.event.emit(event)
            except Exception as exc: self.failed.emit(str(exc))
            finally: self.finished.emit()

    class UploadWorker(QObject):
        result = Signal(dict); failed = Signal(str); finished = Signal()
        def __init__(self, path: str, data_type: str = "expression") -> None:
            super().__init__(); self.path = path; self.data_type = data_type
        def execute(self) -> None:
            try:
                suffix = Path(self.path).suffix.lower()
                if self.data_type == "promoter" or suffix in {".fasta", ".fa"}:
                    records, current, sequence = {}, None, []
                    for line in Path(self.path).read_text(encoding="utf-8").splitlines():
                        if line.startswith(">"):
                            if current: records[current] = "".join(sequence)
                            current, sequence = line[1:].split()[0], []
                        elif current: sequence.append(line.strip())
                    if current: records[current] = "".join(sequence)
                    self.result.emit({"kind": "promoter", "payload": records, "path": self.path,
                                      "detection_reason": "显式选择或 FASTA 后缀"})
                elif self.data_type == "metadata":
                    from phyto_reason.ingestion.parsers.metadata_parser import parse_metadata
                    parsed_md = parse_metadata(self.path)
                    self.result.emit({
                        "kind": "metadata",
                        "payload": {
                            "source_file": Path(self.path).name,
                            "sample_ids": parsed_md.sample_ids,
                            "columns_present": parsed_md.columns_present,
                            "columns_missing": parsed_md.columns_missing,
                            "samples": [sample.model_dump() for sample in parsed_md.samples],
                        },
                        "path": self.path, "warnings": parsed_md.warnings,
                        "detection_reason": "显式选择分组表",
                    })
                else:
                    # Explicit type is authoritative. The filename heuristic is
                    # only a fallback for the default expression selection.
                    if self.data_type in {"metabolite", "metabolite_long"}:
                        from phyto_reason.ingestion.parsers.metabolite_parser import parse_metabolite
                        parsed = parse_metabolite(self.path)
                        reason = ("显式选择长表（m/z-RT-intensity，已透视为宽表）"
                                  if self.data_type == "metabolite_long" else "显式选择代谢物矩阵")
                        self.result.emit({"kind": "metabolite", "payload": parsed.to_workflow_dict(),
                                          "path": self.path, "provenance": parsed.provenance,
                                          "detection_reason": reason})
                        return
                    name = Path(self.path).name.lower()
                    stem = Path(self.path).stem.lower()
                    metabolite_hint = (
                        any(token in name for token in ("metabolite", "metabolom"))
                        or stem in {"meta", "metab"}
                    )
                    if metabolite_hint:
                        try:
                            from phyto_reason.ingestion.parsers.metabolite_parser import parse_metabolite
                            parsed = parse_metabolite(self.path)
                            self.result.emit({"kind": "metabolite", "payload": parsed.to_workflow_dict(),
                                              "path": self.path, "provenance": parsed.provenance,
                                              "detection_reason": "表达矩阵解析的文件名 fallback：检测到代谢物关键词"})
                            return
                        except Exception:
                            pass
                    from phyto_reason.ingestion.parsers.expression_parser import parse_expression
                    parsed = parse_expression(self.path)
                    self.result.emit({"kind": "expression", "payload": parsed.to_workflow_dict(),
                                      "path": self.path, "provenance": parsed.provenance,
                                      "detection_reason": "显式选择表达矩阵"})
            except Exception as exc: self.failed.emit(str(exc))
            finally: self.finished.emit()

    class Window(QMainWindow):
        nav_items = ["总览", "分析", "假设池", "图库", "设置"]
        sources = [("文献", "7", "Europe PMC、PubMed、OpenAlex、bioRxiv、Crossref、arXiv、Semantic Scholar"), ("文献（订阅）", "3", "CAB Abstracts、Scopus、CNKI"), ("通用植物", "4", "PlantTFDB、JASPAR、KEGG、Ensembl Plants"), ("药用植物", "1", "CMAUP v2.0 本地数据库"), ("基因组", "1", "NCBI Datasets API"), ("序列比对", "1", "NCBI remote BLAST"), ("质谱", "2", "MassBank、MoNA"), ("化合物", "1", "PubChem")]
        def __init__(self) -> None:
            super().__init__(); from phyto_reason.desktop.workbench_db import WorkbenchDB; self.setWindowTitle("PhytoReason"); self.resize(1380, 860); self.db = WorkbenchDB(); self.agent = None; self.thread = self.worker = self.upload_thread = self.upload_worker = None; self.session_id = ""; self.response_buffer = ""; self.page_index = {name: i for i, name in enumerate(self.nav_items)}; self._build()
        def text(self, value: str, name: str = "") -> QLabel: label = QLabel(value); label.setObjectName(name); return label
        def _build(self) -> None:
            from phyto_reason.config.llm_config import mask_key
            settings = _legacy_llm_settings(self.db)
            self.base_url = QLineEdit(settings["base_url"]); self.model = QLineEdit(settings["model"]); self.key = QLineEdit(); self.key.setEchoMode(QLineEdit.EchoMode.Password); self.key.setPlaceholderText(f"已配置（{mask_key(settings['api_key'])}），留空保留" if settings["api_key"] else "输入 API Key")
            root = QWidget(); outer = QVBoxLayout(root); outer.setContentsMargins(0,0,0,0); outer.setSpacing(0); top = QFrame(); top.setObjectName("topbar"); tl = QHBoxLayout(top); tl.setContentsMargins(28,14,28,14); tl.addWidget(self.text("PhytoReason", "logo")); tl.addWidget(self.text("植物科研分析工作台", "subtle")); tl.addStretch(); tl.addWidget(self.text("项目", "subtle")); self.project_choice = QComboBox(); self.project_choice.setMinimumWidth(190); tl.addWidget(self.project_choice); self.search_box = QLineEdit(); self.search_box.setPlaceholderText("搜索项目、假设或文献..."); self.search_box.setMinimumWidth(220); tl.addWidget(self.search_box); self.model_choice = QComboBox(); self.model_choice.addItems(["DeepSeek", "OpenAI", "自定义"]); self.model_choice.currentTextChanged.connect(self._preset_changed); tl.addWidget(self.model_choice); outer.addWidget(top)
            body = QHBoxLayout(); body.setContentsMargins(0,0,0,0); body.setSpacing(0); outer.addLayout(body,1); nav = QFrame(); nav.setObjectName("nav"); nav.setFixedWidth(210); nl = QVBoxLayout(nav); nl.setContentsMargins(14,22,14,16); nl.addWidget(self.text("工作台", "navTitle")); self.nav_buttons=[]
            for item in self.nav_items:
                button = QPushButton(item); button.setCheckable(True); button.setObjectName("navButton"); button.clicked.connect(lambda _=False, value=item: self._navigate(value)); nl.addWidget(button); self.nav_buttons.append(button)
            nl.addStretch(); nl.addWidget(self.text("本地模式 · 数据保存在本机", "navHint")); body.addWidget(nav); self.pages = QStackedWidget(); self.pages.addWidget(self._dashboard()); self.pages.addWidget(self._analysis()); self.pages.addWidget(self._hypothesis_pool()); self.pages.addWidget(self._gallery()); self.pages.addWidget(self._settings_page()); body.addWidget(self.pages,1); status = QFrame(); status.setObjectName("statusbar"); sl = QHBoxLayout(status); sl.setContentsMargins(24,4,24,4); self.status = self.text("● Agent 就绪", "statusText"); sl.addWidget(self.status); sl.addStretch(); sl.addWidget(self.text("WorkbenchDB · SQLite WAL", "statusText")); outer.addWidget(status); self.setCentralWidget(root); self._refresh_projects(); self._navigate("总览")
            self.setStyleSheet("""QMainWindow,QWidget{background:#f4f3ef;color:#2c2a29;font-family:'Microsoft YaHei UI','PingFang SC';font-size:13px}#topbar{background:#fff;border-bottom:1px solid #e6e3dc}#logo{font-size:19px;font-weight:700}#subtle,.muted{color:rgba(44,42,41,.62)}#nav{background:#fff;border-right:1px solid #e6e3dc}#navTitle{color:rgba(44,42,41,.4);font-size:12px;padding:4px 12px 10px}#navButton{text-align:left;background:transparent;border:0;border-radius:11px;color:rgba(44,42,41,.62);padding:11px 13px}#navButton:hover,#navButton:checked{background:#d4e0ce;color:#4a6353;font-weight:600}#navHint{color:rgba(44,42,41,.4);font-size:11px;padding:8px 10px}#statusbar{height:40px;background:#fff;border-top:1px solid #e6e3dc}#statusText{color:#4a6353;font-size:11px}.card{background:#fff;border:1px solid #e6e3dc;border-radius:16px}.pageTitle{font-size:26px;font-weight:700}.sectionTitle{font-size:14px;font-weight:600}.metric{font-size:32px;font-weight:700}.metricLabel{color:rgba(44,42,41,.62);font-size:12px}#pairTitle{font-size:14.5px;font-weight:600}#ratingBadge{padding:2px 10px;border-radius:7px;font-size:11px;font-weight:600}QLabel[rating="plausible"]{background:#d4e0ce;color:#4a6353}QLabel[rating="weak"]{background:#e6cdb1;color:#8a6437}QLabel[rating="insufficient"]{background:#e6e3dc;color:#75716a}#chip{font-size:10.5px;padding:2px 8px;border-radius:6px;background:#f4f3ef;color:rgba(44,42,41,.62);border:1px solid #e6e3dc}QFrame#darkCard{background:#2c2a29;border:0;border-radius:16px}QFrame#darkCard QLabel{color:#f4f3ef}QPushButton{background:#e6cdb1;color:#2c2a29;border:0;border-radius:11px;padding:9px 16px;font-weight:600}QPushButton:hover{background:#dfc19e}QPushButton#secondary{background:#fff;color:#2c2a29;border:1px solid #e6e3dc}QPushButton#secondary:hover{background:#f4f3ef}QComboBox,QLineEdit,QTextBrowser,QListWidget{background:#fff;border:1px solid #e6e3dc;border-radius:10px;padding:8px}QListWidget::item{padding:9px;border-radius:10px}QListWidget::item:selected{background:#d4e0ce;color:#4a6353}QScrollArea{background:transparent;border:0}""")
            self.setStyleSheet(self.styleSheet() + "#taskDot{color:#e6cdb1;font-size:9px}#taskStatus{color:rgba(244,243,239,.5);font-size:11px}")
        def _card(self,title,value,note):
            f=QFrame(); f.setObjectName("card"); l=QVBoxLayout(f); l.setContentsMargins(18,16,18,16); l.addWidget(self.text(title,"metricLabel")); l.addWidget(self.text(value,"metric")); l.addWidget(self.text(note,"muted")); return f
        def _dashboard(self):
            page=QWidget(); l=QVBoxLayout(page); l.setContentsMargins(32,26,32,24); l.setSpacing(16); l.addWidget(self.text("研究工作台","pageTitle")); l.addWidget(self.text("项目、数据集、候选机制和分析结果总览。","muted")); row=QHBoxLayout(); self.dashboard_cards={}
            for title,key,note in [("研究项目","projects","本地项目"),("数据集","datasets","已登记数据"),("待验证假设","hypotheses","三档评级记录"),("可视化结果","charts","历史图表")]: self.dashboard_cards[key]=self._card(title,"0",note); row.addWidget(self.dashboard_cards[key])
            l.addLayout(row); mid=QHBoxLayout(); trend=ClickableCard(); tv=QVBoxLayout(trend); tv.setContentsMargins(20,18,20,18); tv.addWidget(self.text("假设评级分布","sectionTitle")); self.rating_summary=self.text("","muted"); tv.addWidget(self.rating_summary); self.dist_bar_host=QWidget(); self.dist_bar=QHBoxLayout(self.dist_bar_host); tv.addWidget(self.dist_bar_host); self.rating_rows=QVBoxLayout(); tv.addLayout(self.rating_rows); self.rating_activity=self.text("","muted"); tv.addWidget(self.rating_activity); trend.clicked.connect(lambda:self._navigate("假设池")); mid.addWidget(trend,2); tasks=QFrame(); tasks.setObjectName("darkCard"); q=QVBoxLayout(tasks); q.setContentsMargins(20,18,20,18); q.addWidget(self.text("运行中任务","sectionTitle")); self.task_rows=QVBoxLayout(); q.addLayout(self.task_rows); q.addStretch(); mid.addWidget(tasks,1); l.addLayout(mid); l.addWidget(self.text("最近假设","sectionTitle")); self.recent_hyp_row=QHBoxLayout(); l.addLayout(self.recent_hyp_row); l.addWidget(self.text("最近项目","sectionTitle")); self.project_list=QListWidget(); l.addWidget(self.project_list); self.refresh_dashboard(); return page
        def _analysis(self):
            page=QWidget(); l=QVBoxLayout(page); l.setContentsMargins(24,22,24,20); title=QHBoxLayout(); title.addWidget(self.text("分析","pageTitle")); title.addStretch(); self.analysis_context=self.text("未选择会话","muted"); title.addWidget(self.analysis_context); l.addLayout(title); split=QSplitter(); left=QWidget(); ll=QVBoxLayout(left); ll.addWidget(self.text("会话","sectionTitle")); new=QPushButton("＋ 新建会话"); new.clicked.connect(self._new_session); ll.addWidget(new); self.sessions=QListWidget(); self.sessions.itemClicked.connect(self._session_clicked); ll.addWidget(self.sessions,1); split.addWidget(left); center=QWidget(); cl=QVBoxLayout(center); self.chat=QTextBrowser(); self.chat.setHtml("<h3>开始一次科研分析</h3><p>上传数据或提出植物代谢调控问题。</p>"); cl.addWidget(self.chat,1); controls=QHBoxLayout(); upload=QPushButton("附件"); upload.setObjectName("secondary"); upload.clicked.connect(self._upload); controls.addWidget(upload); self.input=QLineEdit(); self.input.setPlaceholderText("输入科研问题或分析要求..."); self.input.returnPressed.connect(self._send); controls.addWidget(self.input,1); send=QPushButton("发送"); send.clicked.connect(self._send); controls.addWidget(send); cl.addLayout(controls); split.addWidget(center); l.addWidget(split,1); self.evidence_dock=QDockWidget("证据链",self); self.evidence_dock.setAllowedAreas(Qt.DockWidgetArea.RightDockWidgetArea); self.evidence_dock.setWidget(self._evidence_widget()); self.evidence_dock.setFeatures(QDockWidget.DockWidgetFeature.DockWidgetClosable); self.addDockWidget(Qt.DockWidgetArea.RightDockWidgetArea,self.evidence_dock); self.evidence_dock.hide(); return page
        def _evidence_widget(self):
            f=QWidget(); l=QVBoxLayout(f); l.setContentsMargins(14,14,14,14); l.addWidget(self.text("证据链抽屉","sectionTitle")); l.addWidget(self.text("工具名、查询参数、原始返回和来源将在此展开。","muted")); l.addStretch(); return f
        def _hypothesis_pool(self):
            page=QWidget(); outer=QVBoxLayout(page); outer.setContentsMargins(32,26,32,24); outer.addWidget(self.text("假设池","pageTitle")); outer.addWidget(self.text("候选调控机制、置信度和证据徽章。","muted")); filters=QHBoxLayout(); self.rating_filter=QComboBox(); self.rating_filter.addItems(["评级：全部","Plausible","Weak","Insufficient"]); self.rating_filter.currentTextChanged.connect(self.refresh_hypotheses); filters.addWidget(self.rating_filter); self.species_filter=QComboBox(); self.species_filter.addItem("物种：全部"); filters.addWidget(self.species_filter); self.hyp_search=QLineEdit(); self.hyp_search.setPlaceholderText("搜索 TF、代谢物或关键词"); self.hyp_search.textChanged.connect(self.refresh_hypotheses); filters.addWidget(self.hyp_search,1); outer.addLayout(filters); self.hyp_summary=self.text("","muted"); outer.addWidget(self.hyp_summary); scroll=QScrollArea(); scroll.setWidgetResizable(True); body=QWidget(); self.hyp_cards=QVBoxLayout(body); scroll.setWidget(body); outer.addWidget(scroll,1); return page
        def _gallery(self):
            page=QWidget(); l=QVBoxLayout(page); l.setContentsMargins(32,26,32,24); l.addWidget(self.text("图库","pageTitle")); l.addWidget(self.text("历史生成图表按会话和类型筛选。","muted")); empty=QFrame(); empty.setObjectName("card"); b=QVBoxLayout(empty); b.setContentsMargins(24,48,24,48); b.addWidget(self.text("暂无图表，去分析页生成","empty"),alignment=Qt.AlignmentFlag.AlignCenter); go=QPushButton("前往分析"); go.setObjectName("secondary"); go.clicked.connect(lambda:self._navigate("分析")); b.addWidget(go,alignment=Qt.AlignmentFlag.AlignCenter); l.addWidget(empty); l.addStretch(); return page
        def _settings_page(self):
            page=QWidget(); scroll=QScrollArea(); scroll.setWidgetResizable(True); body=QWidget(); l=QVBoxLayout(body); l.setContentsMargins(32,26,32,32); l.addWidget(self.text("设置","pageTitle")); l.addWidget(self.text("API 配置、数据源管理、外观和缓存。","muted")); api=QFrame(); api.setObjectName("card"); form=QFormLayout(api); form.setContentsMargins(20,18,20,18); form.addRow("Base URL",self.base_url); form.addRow("模型",self.model); form.addRow("API Key",self.key); save=QPushButton("保存配置"); save.clicked.connect(self._configure); form.addRow("",save); l.addWidget(api); l.addWidget(self.text("数据源管理","sectionTitle")); grid=QGridLayout(); grid.setSpacing(12); grid_positions={0:(0,0),1:(0,1),2:(1,0),3:(2,0,1,2),4:(3,0),5:(3,1),6:(4,0),7:(4,1)}
            for i,(name,count,detail) in enumerate(self.sources):
                f=QFrame(); f.setObjectName("card"); b=QVBoxLayout(f); b.setContentsMargins(16,14,16,14); b.addWidget(self.text(name,"sectionTitle")); ready=name=="药用植物" and (_PROJECT_ROOT / "cmaup_data").exists(); summary="0 物种 / 0 成分 / 0 靶点" if name=="药用植物" else f"{count} 个来源"; row=QHBoxLayout(); row.addWidget(self.text(f"● {summary} · {'本地数据已就绪' if ready else '待接入'}","muted")); row.addStretch(); toggle=QCheckBox("启用"); toggle.setChecked(ready); row.addWidget(toggle)
                if name=="药用植物": detail_btn=QPushButton("详情"); detail_btn.setObjectName("secondary"); detail_btn.clicked.connect(self._show_cmaup_detail); row.addWidget(detail_btn)
                b.addLayout(row); b.addWidget(self.text(detail,"muted")); grid.addWidget(f,*grid_positions[i])
            l.addLayout(grid); cache=QFrame(); cache.setObjectName("card"); cf=QFormLayout(cache); cf.setContentsMargins(20,18,20,18); theme=QComboBox(); theme.addItem("暖米自然（当前）"); cf.addRow("外观",theme); self.cache_ttl=QComboBox(); self.cache_ttl.addItems(["24 小时","7 天","30 天"]); self.cache_ttl.setCurrentText(self.db.get_setting("cache_ttl","24 小时")); self.cache_ttl.currentTextChanged.connect(lambda value:self.db.set_setting("cache_ttl",value)); cf.addRow("缓存时长",self.cache_ttl); from phyto_reason.platform_paths import cache_dir as _cache_dir; cf.addRow("缓存目录",self.text(str(_cache_dir()),"muted")); l.addWidget(cache); l.addStretch(); scroll.setWidget(body); wrapper=QVBoxLayout(page); wrapper.setContentsMargins(0,0,0,0); wrapper.addWidget(scroll); return page
        def _show_cmaup_detail(self): QMessageBox.information(self,"CMAUP 数据源详情","数据来源：CMAUP v2.0（官方 8 个数据文件）\n许可：以 Download_Readme 所载为准\n本地状态：原始文件未改动；规范化副本位于 cmaup_data/derived/")
        def _navigate(self,text):
            self.pages.setCurrentIndex(self.page_index[text]); [b.setChecked(i==self.page_index[text]) for i,b in enumerate(self.nav_buttons)]; self.status.setText(f"● Agent 就绪 · 当前页面：{text}");
            if text=="总览": self.refresh_dashboard()
            elif text=="假设池": self.refresh_hypotheses()
            elif text=="分析": self._refresh_sessions()
            if text not in ("分析","假设池"): self.evidence_dock.hide()
        def _refresh_projects(self):
            self.project_choice.clear(); self.project_choice.addItem("未选择项目",""); [self.project_choice.addItem(row["name"],row["project_id"]) for row in self.db.list_projects()]
        def refresh_dashboard(self):
            if not getattr(self, "_project_bound", False):
                self.project_list.itemDoubleClicked.connect(lambda _item: self._navigate("分析"))
                self._project_bound = True
            counts=self.db.counts(); [widget.findChildren(QLabel,"metric")[0].setText(str(counts[key])) for key,widget in self.dashboard_cards.items()]; dist=self.db.rating_distribution(); total=sum(dist.values()); self.rating_summary.setText(" · ".join(f"{k} {dist[k]}" for k in ("Plausible","Weak","Insufficient")))
            while self.dist_bar.count(): item=self.dist_bar.takeAt(0); item.widget() and item.widget().deleteLater()
            self.dist_bar_host.setVisible(total>0)
            for key,color in [("Plausible","#d4e0ce"),("Weak","#e6cdb1"),("Insufficient","#e6e3dc")]:
                if dist[key]: seg=QFrame(); seg.setFixedHeight(14); seg.setStyleSheet(f"background:{color};border-radius:7px;"); self.dist_bar.addWidget(seg,dist[key])
            while self.rating_rows.count(): item=self.rating_rows.takeAt(0); item.widget() and item.widget().deleteLater()
            for key in ("Plausible","Weak","Insufficient"):
                host=QWidget(); line=QHBoxLayout(host); line.setContentsMargins(0,0,0,0); badge=QLabel(key); badge.setObjectName("ratingBadge"); badge.setProperty("rating",HypothesisCard.RATING_CLASS[key]); badge.setAlignment(Qt.AlignmentFlag.AlignCenter); badge.setFixedWidth(100); line.addWidget(badge); line.addWidget(self.text(f"{dist[key]} 个","muted")); line.addStretch(); self.rating_rows.addWidget(host)
            recent=self.db.created_since((datetime.now(timezone.utc)-timedelta(days=30)).isoformat()); self.rating_activity.setText(f"累计 {total} 个假设 · 近 30 天新增 {recent}" if total else "暂无假设——在分析页生成第一个假设后，这里显示评级分布")
            while self.task_rows.count(): item=self.task_rows.takeAt(0); item.widget() and item.widget().deleteLater()
            tasks=self.db.list_tasks();
            if not tasks: self.task_rows.addWidget(self.text("暂无运行中任务","muted"))
            for row in tasks: host=QWidget(); line=QHBoxLayout(host); line.setContentsMargins(0,0,0,0); line.addWidget(self.text("●","taskDot")); line.addWidget(self.text(str(row["task_type"]),"taskName")); line.addStretch(); line.addWidget(self.text(str(row["status"]),"taskStatus")); self.task_rows.addWidget(host)
            while self.recent_hyp_row.count(): item=self.recent_hyp_row.takeAt(0); item.widget() and item.widget().deleteLater()
            rows=self.db.list_hypotheses()[:5];
            if not rows: self.recent_hyp_row.addWidget(self.text("暂无假设——在分析页生成第一个假设","muted"))
            for row in rows: card=HypothesisCard(row["hypothesis_id"],row["statement"],row["rating"],row["confidence"]); card.setMaximumWidth(320); card.card_clicked.connect(lambda _id:self._navigate("分析")); self.recent_hyp_row.addWidget(card)
            self.project_list.clear(); projects=self.db.list_projects(); [self.project_list.addItem(row["name"]) for row in projects] or self.project_list.addItem("暂无项目——在分析页开始你的第一个分析")
        def refresh_hypotheses(self):
            while self.hyp_cards.count(): item=self.hyp_cards.takeAt(0); item.widget() and item.widget().deleteLater()
            rows=self.db.list_hypotheses(); rating=self.rating_filter.currentText(); query=self.hyp_search.text().lower(); dist=self.db.rating_distribution(); self.hyp_summary.setText(f"共 {len(rows)} 个候选假设：Plausible {dist['Plausible']} · Weak {dist['Weak']} · Insufficient {dist['Insufficient']}"); rows=[r for r in rows if (rating=="评级：全部" or r["rating"]==rating) and (not query or query in r["statement"].lower())]
            if not rows: self.hyp_cards.addWidget(self.text("暂无假设，在分析页生成第一个假设","muted"))
            for row in rows: card=HypothesisCard(row["hypothesis_id"],row["statement"],row["rating"],row["confidence"],_parse_chips(row["evidence_json"]),True); card.card_clicked.connect(self._show_evidence); self.hyp_cards.addWidget(card)
            self.hyp_cards.addStretch()
        def _show_evidence(self,hypothesis_id=""):
            if hypothesis_id: self.status.setText(f"证据链：{hypothesis_id}")
            self.evidence_dock.show(); self.evidence_dock.raise_()
        def _refresh_sessions(self): self.sessions.clear(); rows=self.db.list_sessions(); [self.sessions.addItem(row["title"] or row["session_id"]) for row in rows] or self.sessions.addItem("新建分析"); self.analysis_context.setText(self.session_id or "未选择会话")
        def _new_session(self): self.session_id=self.db.create_session("新建分析",self.project_choice.currentData() or None); self._refresh_sessions(); self.chat.clear(); self.chat.append("<h3>新建分析会话</h3><p>请输入问题或上传数据。</p>")
        def _session_clicked(self,item): self.analysis_context.setText(item.text())
        def _preset_changed(self,name):
            presets={"DeepSeek":("https://api.deepseek.com/v1","deepseek-chat"),"OpenAI":("https://api.openai.com/v1","gpt-4o-mini")}; name in presets and (self.base_url.setText(presets[name][0]),self.model.setText(presets[name][1]))
        def _configure(self):
            from phyto_reason.config.llm_config import save_config
            provider = {"DeepSeek": "deepseek", "OpenAI": "openai"}.get(self.model_choice.currentText(), "custom")
            payload = {"provider": provider, "base_url": self.base_url.text().strip(), "model": self.model.text().strip(), "api_key": self.key.text()}
            try:
                save_config(payload)
            except ValueError as exc:
                QMessageBox.warning(self, "配置无效", str(exc)); return
            self.db.set_setting("base_url", self.base_url.text().strip()); self.db.set_setting("model", self.model.text().strip()); self.db.set_setting("api_key", self.key.text()); self.key.setText(""); from phyto_reason.config.llm_config import mask_key; self.key.setPlaceholderText("已配置，留空保留"); from phyto_reason.agent.orchestrator import AgentOrchestrator; self.agent=AgentOrchestrator(); self.status.setText("● Agent 就绪 · 模型配置已保存（config.json）")
        def _send(self):
            question=self.input.text().strip()
            if not question or self.thread is not None: return
            if self.agent is None: self._configure()
            if not self.session_id: self._new_session()
            self.chat.append(f"<p><b>你：</b>{html.escape(question)}</p>"); self.response_buffer=""; self.input.clear(); self.status.setText("正在分析..."); self.thread=QThread(); self.worker=AgentWorker(self.agent,question,self.session_id); self.worker.moveToThread(self.thread); self.thread.started.connect(self.worker.execute); self.worker.event.connect(self._event); self.worker.failed.connect(lambda error:self.chat.append(f"<p style='color:#b86b4b'>错误：{html.escape(error)}</p>")); self.worker.finished.connect(self.thread.quit); self.thread.finished.connect(self._clear); self.thread.start()
        def _event(self,event):
            kind=event.get("type")
            if kind=="chunk": self.response_buffer+=event.get("content","")
            elif kind=="tool_start": self.chat.append(f"<span style='background:#f4f3ef;border:1px solid #e6e3dc;border-radius:20px;padding:4px 12px'>🔧 {html.escape(event.get('tool',''))} ⟳</span>")
            elif kind=="tool_end": self.chat.append(f"<span style='background:#d4e0ce;color:#4a6353;border-radius:20px;padding:4px 12px'>🔧 {html.escape(event.get('tool',''))} ✓</span>")
            elif kind=="done": self.session_id=event.get("session_id",self.session_id); clean=self.response_buffer.replace("**","").replace("---",""); self.chat.append(f"<p><b>PhytoReason：</b></p><div>{html.escape(clean).replace(chr(10),'<br>')}</div>"); self.status.setText("● Agent 就绪")
        def _clear(self): self.thread=self.worker=None; self._refresh_sessions()
        def _upload(self):
            path,_=QFileDialog.getOpenFileName(self,"选择组学文件","","Data files (*.csv *.tsv *.xlsx *.fasta *.fa);;All files (*)")
            if not path or self.upload_thread is not None: return
            labels = ["表达矩阵", "代谢物矩阵", "长表 (m/z-RT-intensity)", "分组表", "启动子序列"]
            label, accepted = QInputDialog.getItem(self, "选择数据类型", "上传内容", labels, 0, False)
            if not accepted:
                return
            data_type = {"表达矩阵": "expression", "代谢物矩阵": "metabolite",
                         "长表 (m/z-RT-intensity)": "metabolite_long",
                         "分组表": "metadata", "启动子序列": "promoter"}[label]
            if not self.session_id: self._new_session()
            self.status.setText("正在读取数据..."); self.upload_thread=QThread(); self.upload_worker=UploadWorker(path, data_type); self.upload_worker.moveToThread(self.upload_thread); self.upload_thread.started.connect(self.upload_worker.execute); self.upload_worker.result.connect(self._upload_done); self.upload_worker.failed.connect(lambda error:self.chat.append(f"<p style='color:#b86b4b'>上传失败：{html.escape(error)}</p>")); self.upload_worker.finished.connect(self.upload_thread.quit); self.upload_thread.finished.connect(self._upload_clear); self.upload_thread.start()
        def _upload_done(self,data):
            if self.agent is None: self._configure()
            if data["kind"]=="metadata":
                self.agent.upload_data(self.session_id, sample_metadata=data["payload"])
            else:
                self.agent.upload_data(self.session_id,
                                       promoter=data["payload"] if data["kind"]=="promoter" else None,
                                       expression=data["payload"] if data["kind"]=="expression" else None,
                                       metabolite=data["payload"] if data["kind"]=="metabolite" else None)
            path=data["path"]; self.db.add_dataset(Path(path).name,data["kind"],path,session_id=self.session_id,project_id=self.project_choice.currentData() or None); reason = data.get("detection_reason", ""); self.chat.append(f"<p style='color:#4a6353'>已上传：{html.escape(Path(path).name)}<br><span style='color:#777'>识别依据：{html.escape(reason)}</span></p>"); self.status.setText("● Agent 就绪")
        def _upload_clear(self): self.upload_thread=self.upload_worker=None

else:
    Window = None


def run_legacy() -> None:
    if QApplication is None: raise SystemExit("Install the desktop extra: pip install 'phyto-reason-agent[desktop]'")
    app=QApplication.instance() or QApplication([]); app.setFont(QFont("Microsoft YaHei UI",10)); window=Window(); window.show(); app.exec()


def run() -> None:
    if "--legacy" in sys.argv[1:]:
        run_legacy()
        return
    from phyto_reason.desktop.webview_host import run as run_webview
    run_webview()


if __name__ == "__main__": run()
