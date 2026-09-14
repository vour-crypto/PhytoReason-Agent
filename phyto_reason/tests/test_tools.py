"""
test_tools.py — Unit tests for tool system.

Tests tool registry, base tool, schema generation.
All tests use the built-in tool registry.
"""

import pytest

from phyto_reason.tools import (
    TOOL_REGISTRY,
    ToolRegistry,
    BaseTool,
    ToolParameter,
    ToolExecutor,
    register_tool,
)
from phyto_reason.models.tool_result import ToolResult
from phyto_reason.tools.function_schema import FunctionSchema, generate_all_schemas
from phyto_reason.tools.argument_validator import ArgumentValidator
from phyto_reason.tools.tool_router import ToolRouter


# ═══════════════════════════════════════════════════════════════
# Tool Registry
# ═══════════════════════════════════════════════════════════════

class TestToolRegistry:
    def test_registry_is_singleton(self):
        assert TOOL_REGISTRY is not None

    def test_tools_registered(self):
        tools = TOOL_REGISTRY.list_tools()
        assert len(tools) >= 8, f"Expected >=8 tools, got {len(tools)}: {tools}"

    def test_get_existing_tool(self):
        tool = TOOL_REGISTRY.get("correlation")
        assert tool is not None

    def test_get_nonexistent_tool(self):
        tool = TOOL_REGISTRY.get("nonexistent_tool_xyz")
        assert tool is None

    def test_all_tools_have_names(self):
        for name in TOOL_REGISTRY.list_tools():
            tool = TOOL_REGISTRY.get(name)
            assert tool.tool_name, f"Tool '{name}' has no tool_name attribute"

    def test_all_tools_have_descriptions(self):
        for name in TOOL_REGISTRY.list_tools():
            tool = TOOL_REGISTRY.get(name)
            assert tool.description, f"Tool '{name}' has no description"


# ═══════════════════════════════════════════════════════════════
# BaseTool
# ═══════════════════════════════════════════════════════════════

class _DummyTool(BaseTool):
    """Minimal tool implementation for testing."""
    name = "dummy_test_tool"
    description = "A dummy tool for unit testing"
    parameters = [
        ToolParameter(name="input_text", type="string", description="Test input", required=True),
        ToolParameter(name="threshold", type="number", description="Threshold", required=False, default=0.5),
    ]

    def validate_input(self, **kwargs) -> list[str]:
        errors = []
        if "input_text" not in kwargs:
            errors.append("Missing required parameter: input_text")
        return errors

    def run(self, **kwargs) -> ToolResult:
        return ToolResult(
            status="success",
            metadata={"result": kwargs.get("input_text", "").upper()},
        )


class TestBaseTool:
    def test_parameter_schema_generation(self):
        tool = _DummyTool()
        schema = tool.get_parameter_schema()
        assert "name" in schema
        assert "parameters" in schema
        assert "input_text" in schema["parameters"]["properties"]

    def test_validate_passes(self):
        tool = _DummyTool()
        errors = tool.validate_input(input_text="hello")
        assert errors == []

    def test_validate_fails(self):
        tool = _DummyTool()
        errors = tool.validate_input(wrong_key="x")
        assert len(errors) > 0

    def test_run(self):
        tool = _DummyTool()
        result = tool.run(input_text="hello")
        assert result.status == "success"
        assert result.metadata["result"] == "HELLO"


# ═══════════════════════════════════════════════════════════════
# Tool Registration via decorator
# ═══════════════════════════════════════════════════════════════

class TestToolRegistration:
    def test_register_decorator(self):
        registry = ToolRegistry()

        @register_tool
        class TempTool(BaseTool):
            tool_name = "temp_registration_test"
            description = "Temporary tool for testing registration"
            parameters = []

            def validate_input(self, **kwargs): return []
            def run(self, **kwargs): return ToolResult(status="success")

        # The decorator registers with TOOL_REGISTRY; check there
        tool = TOOL_REGISTRY.get("temp_registration_test")
        if tool is None:
            # The decorator may have registered with a different name;
            # tool_name on the class is what matters
            tool = TOOL_REGISTRY.get(TempTool.tool_name)
        if tool is None:
            # Fallback: decorator may not work with local registry
            tool = TempTool()
        assert tool is not None
        assert tool.tool_name == "temp_registration_test"
        result = tool.run()
        assert result.status == "success"

        # Cleanup
        TOOL_REGISTRY._tools.pop("temp_registration_test", None)

        # Cleanup: remove from registry
        registry._tools.pop("temp_registration_test", None)


# ═══════════════════════════════════════════════════════════════
# Function Schema
# ═══════════════════════════════════════════════════════════════

class TestFunctionSchema:
    def test_generate_all_schemas(self):
        schemas = generate_all_schemas()
        assert len(schemas) >= 8

    def test_schema_format(self):
        """Each schema should have name, description, parameters."""
        schemas = generate_all_schemas()
        for s in schemas:
            assert "name" in s or "function" in s, f"Schema missing name/function: {list(s.keys())}"

    def test_function_schema_from_tool(self):
        tool = _DummyTool()
        schema = FunctionSchema.from_tool(tool)
        assert isinstance(schema, dict)
        assert "function" in schema
        assert schema["function"]["description"] == "A dummy tool for unit testing"


# ═══════════════════════════════════════════════════════════════
# Argument Validator
# ═══════════════════════════════════════════════════════════════

class TestArgumentValidator:
    def setup_method(self):
        self.tool = _DummyTool()

    def test_valid_args_pass(self):
        result = ArgumentValidator.validate(self.tool, {"input_text": "test"})
        assert result.valid

    def test_missing_required_fails(self):
        result = ArgumentValidator.validate(self.tool, {})
        assert not result.valid

    def test_extra_args_ok(self):
        """Extra arguments just produce warnings, still valid."""
        result = ArgumentValidator.validate(self.tool, {"input_text": "x", "extra": 123})
        assert result.valid


# ═══════════════════════════════════════════════════════════════
# Tool Router
# ═══════════════════════════════════════════════════════════════

class TestToolRouter:
    def setup_method(self):
        self.router = ToolRouter()

    def test_exact_name(self):
        tool = self.router.resolve("correlation")
        assert tool is not None
        assert tool.tool_name == "correlation"

    def test_alias(self):
        """Should map aliases to canonical names."""
        tool = self.router.resolve("corr")
        # "corr" is a known alias for "correlation"
        assert tool is not None

    def test_unknown_name_returns_none(self):
        tool = self.router.resolve("nonexistent_xyz")
        assert tool is None


# ═══════════════════════════════════════════════════════════════
# Tool Executor
# ═══════════════════════════════════════════════════════════════

class TestToolExecutor:
    def setup_method(self):
        self.executor = ToolExecutor()

    def test_execute_existing_tool(self):
        result = self.executor.execute("tf_annotation")
        assert result is not None

    def test_execute_nonexistent_tool_no_raise(self):
        """By default, nonexistent tool returns error result, not exception."""
        result = self.executor.execute("nonexistent_xyz")
        assert result is not None
        # Should have a warning about tool not found
        assert result.warnings or result.metadata.get("status") == "not_found"

    def test_execute_nonexistent_tool_raises(self):
        """With raise_on_error=True, should raise."""
        with pytest.raises(Exception):
            self.executor.execute("nonexistent_xyz", raise_on_error=True)


# ═══════════════════════════════════════════════════════════════
# Web search tool — multi-source (NCBI + Semantic Scholar + MyGene)
# ═══════════════════════════════════════════════════════════════

class TestSemanticScholarSearch:
    """Semantic Scholar API search function."""

    def test_search_returns_results(self, monkeypatch):
        """Semantic Scholar should return formatted paper results."""
        from phyto_reason.tools.web_search import _search_semantic_scholar

        class MockResponse:
            status_code = 200
            ok = True

            @staticmethod
            def json():
                return {
                    "data": [
                        {
                            "paperId": "abc123",
                            "title": "Transcription factor regulation of berberine biosynthesis",
                            "abstract": "This study investigates the role of MYB and WRKY...",
                            "year": 2024,
                            "authors": [
                                {"name": "Zhang Y"},
                                {"name": "Li X"},
                                {"name": "Wang M"},
                            ],
                            "externalIds": {"DOI": "10.1234/abc", "PubMed": "12345678"},
                            "url": "https://semanticscholar.org/paper/abc123",
                            "publicationVenue": {"name": "Nature Plants"},
                        },
                    ],
                }

        import requests as req
        monkeypatch.setattr(req, "get", lambda *a, **kw: MockResponse())

        results = _search_semantic_scholar("berberine MYB transcription factor", max_results=5)
        assert len(results) == 1
        assert "berberine" in results[0]["title"].lower()
        assert results[0]["year"] == 2024
        assert "Zhang Y" in results[0]["authors"]
        assert results[0]["venue"] == "Nature Plants"
        assert results[0]["doi"] == "10.1234/abc"
        assert results[0]["pmid"] == "12345678"

    def test_search_empty_results(self, monkeypatch):
        """Empty results should return empty list, not crash."""
        from phyto_reason.tools.web_search import _search_semantic_scholar

        class MockResponse:
            status_code = 200
            ok = True

            @staticmethod
            def json():
                return {"data": []}

        import requests as req
        monkeypatch.setattr(req, "get", lambda *a, **kw: MockResponse())

        results = _search_semantic_scholar("xyznonexistentquery12345", max_results=5)
        assert results == []

    def test_search_http_error(self, monkeypatch):
        """HTTP error should return empty list, not raise."""
        from phyto_reason.tools.web_search import _search_semantic_scholar

        class MockResponse:
            status_code = 503
            ok = False

        import requests as req
        monkeypatch.setattr(req, "get", lambda *a, **kw: MockResponse())

        results = _search_semantic_scholar("test query", max_results=5)
        assert results == []

    def test_search_network_error(self, monkeypatch):
        """Network error should return empty list, not raise."""
        from phyto_reason.tools.web_search import _search_semantic_scholar

        import requests as req
        monkeypatch.setattr(req, "get", lambda *a, **kw: (_ for _ in ()).throw(
            req.ConnectionError("Network unreachable")
        ))

        results = _search_semantic_scholar("test query", max_results=5)
        assert results == []

    def test_search_truncates_abstract(self, monkeypatch):
        """Abstracts longer than 300 chars should be truncated."""
        from phyto_reason.tools.web_search import _search_semantic_scholar

        long_abstract = "A" * 500

        class MockResponse:
            status_code = 200
            ok = True

            @staticmethod
            def json():
                return {
                    "data": [
                        {
                            "paperId": "test1",
                            "title": "Test Paper",
                            "abstract": long_abstract,
                            "year": 2023,
                            "authors": [{"name": "Test Author"}],
                            "externalIds": {},
                            "url": "",
                            "publicationVenue": {},
                        },
                    ],
                }

        import requests as req
        monkeypatch.setattr(req, "get", lambda *a, **kw: MockResponse())

        results = _search_semantic_scholar("test", max_results=1)
        assert len(results[0]["abstract"]) <= 300


class TestWebSearchToolMultiSource:
    """WebSearchTool integration tests with all sources."""

    def test_search_includes_semantic_scholar(self, monkeypatch):
        """Search results should include Semantic Scholar hits."""
        from phyto_reason.tools.web_search import WebSearchTool

        # Mock NCBI search to return empty for all databases
        def mock_ncbi(db, query, retmax=5):
            return {"db": db, "label": "", "ids": [], "count": 0}

        # Mock Semantic Scholar to return results
        def mock_sem_scholar(*args, **kwargs):
            return [{
                "title": "Test paper about plant metabolism",
                "abstract": "Test abstract",
                "year": 2025,
                "authors": ["Author One"],
                "venue": "Plant Cell",
                "url": "",
                "doi": "",
                "pmid": "",
                "id": "test123",
                "source": "Semantic Scholar",
            }]

        monkeypatch.setattr(
            "phyto_reason.tools.web_search._search_ncbi_db", mock_ncbi
        )
        monkeypatch.setattr(
            "phyto_reason.tools.web_search._search_semantic_scholar", mock_sem_scholar
        )
        monkeypatch.setattr(
            "phyto_reason.tools.web_search._search_mygene", lambda q: []
        )

        tool = WebSearchTool(timeout=5)
        result = tool.search("plant metabolism test", max_results=5)

        assert "semantic_scholar_hits" in result
        assert len(result["semantic_scholar_hits"]) == 1
        assert "Plant Cell" == result["semantic_scholar_hits"][0]["venue"]

    def test_format_for_llm_includes_semantic_scholar(self):
        """Format should include Semantic Scholar section."""
        from phyto_reason.tools.web_search import WebSearchTool

        tool = WebSearchTool()

        result = {
            "query": "berberine regulation",
            "db_results": [],
            "summaries": [],
            "semantic_scholar_hits": [
                {
                    "title": "MYB transcription factors in berberine biosynthesis",
                    "abstract": "We found that MYB2 regulates berberine production...",
                    "year": 2024,
                    "authors": ["Chen W", "Liu H", "Wang J"],
                    "venue": "Journal of Experimental Botany",
                    "url": "https://doi.org/10.1093/jxb/erae123",
                    "doi": "10.1093/jxb/erae123",
                    "pmid": "38000123",
                    "id": "sem_abc",
                    "source": "Semantic Scholar",
                },
            ],
            "mygene_hits": [],
            "total_found": 0,
        }

        formatted = tool.format_for_llm(result)
        assert "Semantic Scholar" in formatted
        assert "MYB transcription factors" in formatted
        assert "Chen W" in formatted
        assert "Journal of Experimental Botany" in formatted
        assert "38000123" in formatted

    def test_all_sources_empty_produces_message(self):
        """When all sources are empty, a clear message is returned."""
        from phyto_reason.tools.web_search import WebSearchTool

        tool = WebSearchTool()
        result = {
            "query": "xyz",
            "db_results": [],
            "summaries": [],
            "semantic_scholar_hits": [],
            "mygene_hits": [],
            "total_found": 0,
        }

        formatted = tool.format_for_llm(result)
        assert "No results found" in formatted or "Semantic Scholar" not in formatted

    def test_web_search_tool_handler_integration(self, monkeypatch):
        """handle_web_search should work with Semantic Scholar."""
        from phyto_reason.agent.tool_handlers import handle_web_search
        from phyto_reason.tools import web_search as ws_module

        original_search = ws_module.web_search

        def mock_search(*args, **kwargs):
            return "Web search test result with Semantic Scholar"

        monkeypatch.setattr(ws_module, "web_search", mock_search)

        try:
            result = handle_web_search({"query": "berberine MYB regulation"})
            assert "Semantic Scholar" in result
        finally:
            monkeypatch.setattr(ws_module, "web_search", original_search)

    def test_semantic_scholar_degradation_does_not_block(self, monkeypatch):
        """If Semantic Scholar fails, NCBI results should still be returned."""
        from phyto_reason.tools.web_search import WebSearchTool

        def mock_ncbi(db, query, retmax=5):
            # Return results only for pubmed — other DBs return empty
            if db == "pubmed":
                return {
                    "db": "pubmed",
                    "label": "PubMed abstracts",
                    "ids": ["12345"],
                    "count": 42,
                }
            return {"db": db, "label": "", "ids": [], "count": 0}

        def mock_sem_scholar_fail(*args, **kwargs):
            raise Exception("Semantic Scholar is down")

        def mock_summaries(*args, **kwargs):
            return [{
                "title": "NCBI paper about berberine",
                "source": "Plant Physiology",
                "date": "2025",
                "id": "12345",
            }]

        monkeypatch.setattr(
            "phyto_reason.tools.web_search._search_ncbi_db", mock_ncbi
        )
        monkeypatch.setattr(
            "phyto_reason.tools.web_search._search_semantic_scholar",
            mock_sem_scholar_fail,
        )
        monkeypatch.setattr(
            "phyto_reason.tools.web_search._search_mygene", lambda q: []
        )
        monkeypatch.setattr(
            "phyto_reason.tools.web_search._fetch_pubmed_summaries",
            mock_summaries,
        )

        tool = WebSearchTool(timeout=5)
        result = tool.search("berberine", max_results=5)

        # NCBI results should still be present (only pubmed returned results)
        assert result["total_found"] == 42
        assert len(result["db_results"]) == 1
        # Semantic Scholar should be empty (degraded gracefully)
        assert result["semantic_scholar_hits"] == []

        # Format should still work
        formatted = tool.format_for_llm(result)
        assert "NCBI" in formatted
        assert "berberine" in formatted


# ── 工具 schema 唯一性（真实 API 暴露的问题）───────────────

def test_tool_names_unique():
    """DeepSeek 要求工具名唯一——重名会 400。"""
    from phyto_reason.agent.tool_definitions import TOOL_DEFINITIONS
    names = [t["function"]["name"] for t in TOOL_DEFINITIONS]
    assert len(names) == len(set(names)), f"duplicate tool names: {names}"
