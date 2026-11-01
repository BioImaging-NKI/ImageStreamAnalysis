"""Command line: imagestream-cli run settings.toml | run-project | overview | gui | example-settings"""

from __future__ import annotations

import argparse
import logging
import sys

from .config import Settings


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="imagestream-cli", description="Analyze ImageStream event stacks")
    sub = parser.add_subparsers(dest="command", required=True)
    run_p = sub.add_parser("run", help="analyze the image given in a settings file")
    run_p.add_argument("settings", help="settings .toml file")
    run_p.add_argument("--max-events", type=int, help="override input.max_events")
    sub.add_parser("example-settings", help="print a settings file with all defaults")
    proj_p = sub.add_parser("run-project", help="analyze all datasets of a project (project.toml or its folder)")
    proj_p.add_argument("project", help="project.toml, or the project folder")
    over_p = sub.add_parser("overview", help="write the overview tables of a project again from its results (for "
                                             "example after changing the conditions), without analysing")
    over_p.add_argument("project", help="project.toml, or the project folder")
    gui_p = sub.add_parser("gui", help="open the import and run page in the browser")
    gui_p.add_argument("--port", type=int, default=0, help="port (default: a free one)")
    gui_p.add_argument("--no-browser", action="store_true", help="do not open the browser")
    args = parser.parse_args(argv)

    if args.command == "example-settings":
        sys.stdout.write(Settings().to_toml())
        return 0

    if args.command == "gui":
        from .gui import serve

        serve(args.port, open_browser=not args.no_browser)
        return 0

    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    if args.command == "run-project":
        from .project import run_project

        run_project(args.project)
        return 0
    if args.command == "overview":
        from pathlib import Path

        from .project import write_overview

        folder = Path(args.project)
        write_overview(folder.parent if folder.is_file() else folder)
        return 0

    from .pipeline import run  # imports scikit-learn: only when needed

    settings = Settings.load(args.settings)
    if args.max_events is not None:
        settings.input.max_events = args.max_events
    run(settings)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
