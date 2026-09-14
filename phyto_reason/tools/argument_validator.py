"""
argument_validator.py — 参数验证层。

基于 tool 的 parameter schema 对传入参数进行类型检查。
"""

from __future__ import annotations

from typing import Any

from pydantic import BaseModel


class ValidationResult(BaseModel):
    valid: bool = True
    errors: list[str] = []
    warnings: list[str] = []


class ArgumentValidator:
    """工具参数验证器。"""

    TYPE_CHECKS = {
        "string": lambda v: isinstance(v, str),
        "number": lambda v: isinstance(v, (int, float)),
        "integer": lambda v: isinstance(v, int),
        "boolean": lambda v: isinstance(v, bool),
        "array": lambda v: isinstance(v, (list, tuple)),
        "object": lambda v: isinstance(v, dict),
    }

    @staticmethod
    def validate(tool, arguments: dict[str, Any]) -> ValidationResult:
        """验证参数是否符合工具的 parameter schema。"""
        errors: list[str] = []
        warnings: list[str] = []
        params = {p.name: p for p in (tool.parameters or [])}

        for p_name, p in params.items():
            if p_name not in arguments or arguments[p_name] is None:
                if p.required:
                    errors.append(f"缺少必需参数: '{p_name}'")
                continue

            value = arguments[p_name]
            check = ArgumentValidator.TYPE_CHECKS.get(p.type)
            if check and not check(value):
                errors.append(
                    f"参数 '{p_name}' 类型错误: 期望 {p.type}, "
                    f"实际 {type(value).__name__}"
                )

        for arg_name in arguments:
            if arg_name not in params:
                warnings.append(f"未知参数: '{arg_name}' (将被忽略)")

        return ValidationResult(
            valid=len(errors) == 0,
            errors=errors,
            warnings=warnings,
        )
