import sys
from collections.abc import Callable, Generator
from typing import Any


def coroutine[F: Callable[..., Generator[Any, Any, Any]]](
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
