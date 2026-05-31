"""
Originally taken from Michael Hall's tbpore https://github.com/mbhall88/tbpore/blob/main/tbpore/external_tools.py

Also used by a variety of other tools (Dnaapler, Plassembler, Pharokka)

"""

import hashlib
import os
import shlex
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import click
from loguru import logger


class ExternalTool:
    def __init__(
        self,
        tool: str,
        input: str,
        output: str,
        params: str,
        logdir: Path,
        env: Optional[Dict[str, str]] = None,
    ):
        logdir = Path(logdir)   # <-- ensure Path
        self.command: List[str] = self._build_command(tool, input, output, params)
        logdir.mkdir(parents=True, exist_ok=True)
        command_hash = hashlib.sha256(self.command_as_str.encode("utf-8")).hexdigest()
        tool_name = Path(tool).name
        logfile_prefix: Path = logdir / f"{tool_name}_{command_hash}"
        self.out_log = f"{logfile_prefix}.out"
        self.err_log = f"{logfile_prefix}.err"
        # Extra env vars merged with os.environ for the subprocess (e.g.
        # CUDA_VISIBLE_DEVICES for multi-GPU foldseek). None == inherit.
        self.env = env

    @property
    def command_as_str(self) -> str:
        return shlex.join(self.command)

    @staticmethod
    def _build_command(tool: str, input: str, output: str, params: str) -> List[str]:
        # note: shlex.join does not allow us to shlex.split() later
        # this is explicitly a " ".join()
        command = " ".join([tool, params, output, input])
        escaped_command = shlex.split(command)
        return escaped_command

    def run(self) -> None:
        with open(self.out_log, "w") as stdout_fh, open(self.err_log, "w") as stderr_fh:
            print(f"Command line: {self.command_as_str}", file=stderr_fh)
            if self.env:
                print(f"Extra env: {self.env}", file=stderr_fh)
            logger.info(f"Started running {self.command_as_str} ...")
            self._run_core(
                self.command,
                stdout_fh=stdout_fh,
                stderr_fh=stderr_fh,
                env=self.env,
            )
            logger.info(f"Done running {self.command_as_str}")

    """
    stream to terminal (aria2c) so the user knows how long it is taking
    """

    def run_stream(self) -> None:
        with open(self.out_log, "w") as stdout_fh, open(self.err_log, "w") as stderr_fh:
            print(f"Command line: {self.command_as_str}", file=stderr_fh)
            logger.info(f"Started running {self.command_as_str} ...")

            process = subprocess.Popen(
                self.command,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                bufsize=1,
                universal_newlines=True,
            )

            for line in process.stdout:
                print(line, end="")  # Live output to terminal
                stdout_fh.write(line)  # Also write to stdout log

            process.stdout.close()
            return_code = process.wait()

            logger.info(f"Done running {self.command_as_str}")

            if return_code != 0:
                raise subprocess.CalledProcessError(return_code, self.command)

    @staticmethod
    def _run_core(
        command: List[str],
        stdout_fh,
        stderr_fh,
        env: Optional[Dict[str, str]] = None,
    ) -> None:
        # Merge `env` (if any) with the inherited environment so PATH,
        # CONDA_PREFIX, etc. survive. None == inherit unchanged.
        merged_env = {**os.environ, **env} if env else None
        subprocess.check_call(
            command, stdout=stdout_fh, stderr=stderr_fh, env=merged_env
        )

    @staticmethod
    def run_tools(
        tools_to_run: Tuple["ExternalTool", ...], ctx: Optional[click.Context] = None
    ) -> None:
        for tool in tools_to_run:
            try:
                tool.run()
            except subprocess.CalledProcessError as error:
                logger.error(
                    f"Error calling {tool.command_as_str} (return code {error.returncode})"
                )
                logger.error(f"Please check stdout log file: {tool.out_log}")
                logger.error(f"Please check stderr log file: {tool.err_log}")
                logger.error("Temporary files are preserved for debugging")
                logger.error("Exiting...")

                if ctx:
                    ctx.exit(1)
                else:
                    sys.exit(1)

    """
    Only one toolf
    """

    @staticmethod
    def run_tool(tool: "ExternalTool", ctx: Optional[click.Context] = None) -> None:
        try:
            tool.run()
        except subprocess.CalledProcessError as error:
            logger.error(
                f"Error calling {tool.command_as_str} (return code {error.returncode})"
            )
            logger.error(f"Please check stdout log file: {tool.out_log}")
            logger.error(f"Please check stderr log file: {tool.err_log}")
            logger.error("Temporary files are preserved for debugging")
            logger.error("Exiting...")

            if ctx:
                ctx.exit(1)
            else:
                sys.exit(1)

    """
    Only download - so can print the aria2c output to screen
    """

    @staticmethod
    def run_download(tool: "ExternalTool", ctx: Optional[click.Context] = None) -> None:
        try:
            tool.run_stream()
        except subprocess.CalledProcessError as error:
            logger.error(
                f"Error calling {tool.command_as_str} (return code {error.returncode})"
            )
            logger.error(f"Please check stdout log file: {tool.out_log}")
            logger.error(f"Please check stderr log file: {tool.err_log}")
            logger.error("Temporary files are preserved for debugging")
            logger.error("Exiting...")

            if ctx:
                ctx.exit(1)
            else:
                sys.exit(1)
