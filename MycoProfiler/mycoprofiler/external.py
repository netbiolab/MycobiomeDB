"""Running external programs, and capturing what was run for the manifest."""

import os
import shutil
import subprocess
import time


class ExternalToolError(Exception):
    """An external command failed. Carries enough context to name the step."""

    def __init__(self, step, command, returncode, log_path=None, stderr_tail=None):
        self.step = step
        self.command = command
        self.returncode = returncode
        self.log_path = log_path
        self.stderr_tail = stderr_tail
        message = ["step '%s' failed: %s exited with status %s" % (step, command[0], returncode)]
        message.append("command: %s" % shell_quote(command))
        if log_path:
            message.append("log: %s" % log_path)
        if stderr_tail:
            message.append("last output:")
            message.extend("  " + line for line in stderr_tail)
        super(ExternalToolError, self).__init__("\n".join(message))


def shell_quote(command):
    """Render an argv list as a copy-pasteable shell command line."""
    try:
        from shlex import quote
    except ImportError:  # pragma: no cover - Python 2 only
        from pipes import quote
    return " ".join(quote(str(arg)) for arg in command)


def find_executable(name, override=None):
    """Locate an executable, honouring an explicit override path."""
    if override:
        if os.path.isfile(override) and os.access(override, os.X_OK):
            return os.path.abspath(override)
        found = shutil.which(override)
        if found:
            return found
        raise ExternalToolError(
            "startup", [override], 127,
            stderr_tail=["'%s' is not an executable file and is not on $PATH" % override],
        )
    found = shutil.which(name)
    if not found:
        raise ExternalToolError(
            "startup", [name], 127,
            stderr_tail=[
                "'%s' was not found on $PATH." % name,
                "Install it (e.g. `conda env create -f environment.yml && conda activate mycoprofiler`)",
                "or pass its path with --kraken2.",
            ],
        )
    return found


def tool_version(executable):
    """Best-effort one-line version string for an external tool."""
    try:
        proc = subprocess.run(
            [executable, "--version"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            universal_newlines=True, timeout=60,
        )
    except (OSError, subprocess.SubprocessError) as exc:
        return "unavailable (%s)" % exc
    text = (proc.stdout or "").strip()
    return text.splitlines()[0] if text else "unavailable (no output)"


def run_command(command, step, log_path, echo=None):
    """Run `command`, tee-ing stdout+stderr to `log_path`.

    Returns the captured output as a list of lines. Raises ExternalToolError on a
    non-zero exit status, so a failed step can never feed the next one.
    """
    if echo:
        echo("$ " + shell_quote(command))
    started = time.time()
    captured = []
    with open(log_path, "w") as log:
        log.write("# step: %s\n" % step)
        log.write("# command: %s\n" % shell_quote(command))
        log.write("# started: %s\n\n" % time.strftime("%Y-%m-%dT%H:%M:%S%z"))
        log.flush()
        # `with Popen(...)` closes the stdout pipe on the way out; iterating it and
        # calling wait() alone leaks the file descriptor.
        with subprocess.Popen(
            command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            universal_newlines=True, bufsize=1,
        ) as proc:
            for line in proc.stdout:
                captured.append(line.rstrip("\n"))
                log.write(line)
            returncode = proc.wait()
        elapsed = time.time() - started
        log.write("\n# exit status: %s\n" % returncode)
        log.write("# elapsed seconds: %.2f\n" % elapsed)

    if returncode != 0:
        raise ExternalToolError(
            step, command, returncode, log_path=log_path, stderr_tail=captured[-15:]
        )
    return captured, elapsed
