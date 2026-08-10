import os

os.chdir(os.path.split(os.path.realpath(__file__))[0])

import argparse
import ast
import platform
import sys
from importlib import metadata
from pathlib import Path

SKIP_DIRS = {
    ".git",
    ".hg",
    ".svn",
    ".idea",
    ".vscode",
    "__pycache__",
    ".pytest_cache",
    ".mypy_cache",
    ".ruff_cache",
    ".venv",
    "venv",
    "env",
    ".env",
    "site-packages",
    "build",
    "dist",
    ".tox",
    ".nox",
}


def get_stdlib_names():
    if hasattr(sys, "stdlib_module_names"):
        return set(sys.stdlib_module_names)
    names = set(sys.builtin_module_names)
    try:
        import pkgutil

        stdlib = Path(os.__file__).resolve().parent
        for m in pkgutil.iter_modules([str(stdlib)]):
            names.add(m.name)
    except Exception:
        pass
    return names


def iter_python_files(root):
    root = Path(root).resolve()
    if root.is_file():
        if root.suffix == ".py":
            yield root
        return
    for p in root.rglob("*.py"):
        if any(part in SKIP_DIRS for part in p.parts):
            continue
        yield p


def top_name(name):
    return name.split(".", 1)[0].strip()


def imports_from_ast(path):
    found = set()
    try:
        text = path.read_text(encoding="utf-8")
    except UnicodeDecodeError:
        text = path.read_text(encoding="utf-8", errors="ignore")
    try:
        tree = ast.parse(text, filename=str(path))
    except SyntaxError as e:
        return found, f"{path}: {e}"

    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for alias in node.names:
                name = top_name(alias.name)
                if name:
                    found.add(name)

        elif isinstance(node, ast.ImportFrom):
            if node.level == 0 and node.module:
                name = top_name(node.module)
                if name:
                    found.add(name)

        elif isinstance(node, ast.Call):
            func = node.func
            is_dynamic_import = (
                isinstance(func, ast.Name) and func.id == "__import__"
            ) or (
                isinstance(func, ast.Attribute)
                and func.attr == "import_module"
                and isinstance(func.value, ast.Name)
                and func.value.id == "importlib"
            )
            if is_dynamic_import and node.args:
                arg = node.args[0]
                if isinstance(arg, ast.Constant) and isinstance(arg.value, str):
                    name = top_name(arg.value)
                    if name:
                        found.add(name)

    return found, None


def discover_local_modules(root):
    root = Path(root).resolve()
    local = set()

    if root.is_file():
        root = root.parent

    for p in root.iterdir():
        if p.name in SKIP_DIRS:
            continue
        if p.is_file() and p.suffix == ".py":
            local.add(p.stem)
        elif p.is_dir():
            if (p / "__init__.py").exists() or any(p.glob("*.py")):
                local.add(p.name)

    src = root / "src"
    if src.is_dir():
        for p in src.iterdir():
            if p.is_file() and p.suffix == ".py":
                local.add(p.stem)
            elif p.is_dir() and ((p / "__init__.py").exists() or any(p.glob("*.py"))):
                local.add(p.name)

    return local


def build_distribution_map():
    try:
        mapping = metadata.packages_distributions()
        return {k: sorted(set(v)) for k, v in mapping.items()}
    except Exception:
        mapping = {}
        for dist in metadata.distributions():
            try:
                names = []
                top_level = dist.read_text("top_level.txt")
                if top_level:
                    names.extend(x.strip() for x in top_level.splitlines() if x.strip())
                dist_name = dist.metadata.get("Name")
                if dist_name:
                    for name in names:
                        mapping.setdefault(name, []).append(dist_name)
            except Exception:
                pass
        return {k: sorted(set(v)) for k, v in mapping.items()}


def normalize_dist_name(name):
    try:
        return metadata.metadata(name).get("Name") or name
    except Exception:
        return name


def resolve_dependencies(imports, dist_map):
    resolved = {}
    unresolved = []
    ambiguous = {}

    for mod in sorted(imports, key=str.lower):
        candidates = dist_map.get(mod, [])
        if not candidates:
            unresolved.append(mod)
            continue

        valid = []
        for dist in candidates:
            try:
                ver = metadata.version(dist)
                valid.append((normalize_dist_name(dist), ver))
            except metadata.PackageNotFoundError:
                pass

        if not valid:
            unresolved.append(mod)
            continue

        unique = {}
        for name, ver in valid:
            unique[name.lower()] = (name, ver)
        valid = sorted(unique.values(), key=lambda x: x[0].lower())

        if len(valid) > 1:
            ambiguous[mod] = valid

        for name, ver in valid:
            resolved[name.lower()] = (name, ver)

    return sorted(resolved.values(), key=lambda x: x[0].lower()), unresolved, ambiguous


def main():
    parser = argparse.ArgumentParser(
        description="Generate requirements.txt from Python imports using versions installed in the current environment."
    )
    parser.add_argument(
        "path",
        nargs="?",
        default=".",
        help="Python file or project directory to scan. Default: current directory.",
    )
    parser.add_argument(
        "-o",
        "--output",
        default="requirements.txt",
        help="Output requirements file. Default: requirements.txt",
    )
    parser.add_argument(
        "--report",
        default="requirements_report.txt",
        help="Diagnostic report file. Default: requirements_report.txt",
    )
    args = parser.parse_args()

    root = Path(args.path).resolve()
    py_files = list(iter_python_files(root))

    imports = set()
    syntax_errors = []

    for path in py_files:
        found, err = imports_from_ast(path)
        imports.update(found)
        if err:
            syntax_errors.append(err)

    stdlib = get_stdlib_names()
    local = discover_local_modules(root)

    third_party_imports = {
        x for x in imports if x not in stdlib and x not in local and x != "__future__"
    }

    dist_map = build_distribution_map()
    requirements, unresolved, ambiguous = resolve_dependencies(
        third_party_imports, dist_map
    )

    output_path = Path(args.output)
    if not output_path.is_absolute():
        output_path = Path.cwd() / output_path

    report_path = Path(args.report)
    if not report_path.is_absolute():
        report_path = Path.cwd() / report_path

    req_lines = [f"{name}=={ver}" for name, ver in requirements]
    output_path.write_text(
        "\n".join(req_lines) + ("\n" if req_lines else ""), encoding="utf-8"
    )

    report = []
    report.append("Import-based dependency report")
    report.append("=" * 32)
    report.append(f"Python: {platform.python_version()}")
    report.append(f"Executable: {sys.executable}")
    report.append(f"Platform: {platform.platform()}")
    report.append(f"Scanned path: {root}")
    report.append(f"Python files: {len(py_files)}")
    report.append("")
    report.append(f"All top-level imports ({len(imports)}):")
    report.extend(f"  {x}" for x in sorted(imports, key=str.lower))
    report.append("")
    report.append(f"Detected local modules ({len(local)}):")
    report.extend(f"  {x}" for x in sorted(local, key=str.lower))
    report.append("")
    report.append(f"Third-party imports ({len(third_party_imports)}):")
    report.extend(f"  {x}" for x in sorted(third_party_imports, key=str.lower))
    report.append("")
    report.append(f"Resolved distributions ({len(requirements)}):")
    report.extend(f"  {name}=={ver}" for name, ver in requirements)
    report.append("")
    report.append(f"Unresolved imports ({len(unresolved)}):")
    report.extend(f"  {x}" for x in unresolved)
    report.append("")
    report.append(f"Ambiguous mappings ({len(ambiguous)}):")
    for mod, values in sorted(ambiguous.items()):
        report.append(f"  {mod}:")
        for name, ver in values:
            report.append(f"    {name}=={ver}")
    report.append("")
    report.append(f"Syntax errors ({len(syntax_errors)}):")
    report.extend(f"  {x}" for x in syntax_errors)

    report_path.write_text("\n".join(report) + "\n", encoding="utf-8")

    print(f"Scanned {len(py_files)} Python file(s).")
    print(f"Found {len(third_party_imports)} third-party import(s).")
    print(f"Wrote {len(requirements)} pinned package(s) to: {output_path}")
    print(f"Report: {report_path}")

    if unresolved:
        print("")
        print("Unresolved imports:")
        for name in unresolved:
            print(f"  - {name}")
        print(
            "These may be local modules, optional dependencies, or packages not installed in this Python environment."
        )


if __name__ == "__main__":
    main()
