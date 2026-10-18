"""Check that the built reference documents every public symbol exactly once."""

import inspect
import sys
import zlib
from pathlib import Path

DOCS_ROOT = Path(__file__).resolve().parent
REPO_ROOT = DOCS_ROOT.parent
DEFAULT_INVENTORY = DOCS_ROOT / "build" / "api" / "objects.inv"


def read_inventory(path: Path) -> list[str]:
    raw = path.read_bytes()
    header = raw.split(b"\n", 4)
    if len(header) < 5:
        raise SystemExit(f"Malformed inventory: {path}")
    lines = zlib.decompress(header[4]).decode("utf-8").splitlines()
    return [line.split(" ", 1)[0] for line in lines if line.strip()]


def public_members(cls: type) -> list[str]:
    """Return the public methods and properties that ``cls`` itself defines."""
    members = []
    for name, member in vars(cls).items():
        if name.startswith("_"):
            continue
        if isinstance(member, classmethod | staticmethod):
            member = member.__func__
        if isinstance(member, property) or inspect.isfunction(member):
            members.append(name)
    return members


def member_report(names: list[str]) -> tuple[list[str], list[str]]:
    from cytetype import CyteType

    missing: list[str] = []
    duplicated: list[str] = []
    for cls in (CyteType,):
        for name in public_members(cls):
            qualified = f"{cls.__name__}.{name}"
            hits = [entry for entry in names if entry.endswith(f".{qualified}")]
            if not hits:
                missing.append(qualified)
            elif len(hits) > 1:
                duplicated.append(f"{qualified}: {', '.join(sorted(hits))}")
    return missing, duplicated


def public_name_report(names: list[str]) -> list[str]:
    import cytetype

    documented = set(names)
    expected = [f"cytetype.{name}" for name in cytetype.__all__]
    return [name for name in expected if name not in documented]


def main() -> None:
    inventory = Path(sys.argv[1]) if len(sys.argv) > 1 else DEFAULT_INVENTORY
    if not inventory.is_file():
        raise SystemExit(f"Build the HTML docs first; no inventory at {inventory}")
    if str(REPO_ROOT) not in sys.path:
        sys.path.insert(0, str(REPO_ROOT))

    names = read_inventory(inventory)
    missing_members, duplicated_members = member_report(names)
    missing_names = public_name_report(names)

    problems: list[str] = []
    if missing_members:
        problems.append(
            "Public CyteType members missing from the reference: "
            + ", ".join(sorted(missing_members))
        )
    if duplicated_members:
        problems.append(
            "Public members documented more than once:\n  "
            + "\n  ".join(sorted(duplicated_members))
        )
    if missing_names:
        problems.append(
            "Public symbols missing from the reference under their public names: "
            + ", ".join(sorted(missing_names))
        )
    if problems:
        raise SystemExit("\n".join(problems))

    print(
        "Reference coverage complete: every cytetype.__all__ symbol is documented under "
        "its public name, and every public CyteType member is documented exactly once"
    )


if __name__ == "__main__":
    main()
