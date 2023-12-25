import logging
from enum import StrEnum, auto
from logging import Formatter
from typing import Dict

from rich.console import Console
from rich.logging import RichHandler

LOG_FORMAT = (
    "%(asctime)s  %(processName)-20s  %(name)-70s  [%(levelname)8s] %(message)s"
)

RICH_CONSOLE = Console(stderr=False)
RICH_HANDLER = RichHandler(
    level=logging.DEBUG,
    log_time_format="[%X]",
    console=RICH_CONSOLE,
    omit_repeated_times=False,
    locals_max_length=None,  # type:ignore
    locals_max_string=None,  # type:ignore
)
RICH_HANDLER.setFormatter(Formatter("%(message)s"))


class LoggingLevel(StrEnum):
    DEBUG = auto()
    INFO = auto()
    WARNING = auto()
    ERROR = auto()
    CRITICAL = auto()


def config_console_debug_logger(name: str):
    log = logging.getLogger(name)
    log.addHandler(RICH_HANDLER)
    log.propagate = False


def log_message(info: Dict, sep: str = " | ") -> str:
    if info is None:
        return  # type:ignore

    repr = []
    for key, value in info.items():
        repr.append(f"[{key}] '{value}'")

    return sep.join(repr)
