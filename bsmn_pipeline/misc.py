import sys
from collections.abc import Callable, Generator
from typing import Any, TypeVar

# TypeVar form (not PEP 695 `def coroutine[F: ...]`) so the module imports on
# Python 3.10 (the cluster's `bp` conda env) as well as 3.12+.
F = TypeVar("F", bound=Callable[..., Generator[Any, Any, Any]])


def coroutine(  # noqa: UP047 - keep TypeVar form for Python 3.10 (cluster bp env)
    func: F,
) -> Callable[..., Generator[Any, Any, Any]]:
    def start(*args: Any, **kwargs: Any) -> Generator[Any, Any, Any]:
        g = func(*args, **kwargs)
        g.__next__()
        return g

    return start


def printer(out: Any) -> None:
    try:
        print(out, flush=True)
    except BrokenPipeError:
        try:
            sys.stdout.close()
        except BrokenPipeError:
            pass
        try:
            sys.stderr.close()
        except BrokenPipeError:
            pass
