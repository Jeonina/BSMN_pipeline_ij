import configparser
import os
import pathlib
import subprocess


def read_config(reference: str = "b37", conda_env: str = "bp") -> configparser.ConfigParser:
    lib_home = os.path.dirname(os.path.realpath(__file__))
    pipe_home = os.path.normpath(lib_home + "/..")
    env_dir = subprocess.check_output(
        f"conda info -e | grep -w ^{conda_env} | awk '{{print $NF}}'",
        shell=True,
        universal_newlines=True,
    ).strip()
    config = configparser.ConfigParser()
    config["PATH"] = {"pipe_home": pipe_home, "env_dir": env_dir}

    config.read(
        pipe_home
        + (
            "/config.hg19.ini"
            if reference == "hg19"
            else "/config.hg38.ini"
            if reference == "hg38"
            else "/config.ini"
        )
    )
    for section in ["TOOLS", "RESOURCES"]:
        for key in config[section]:
            config[section][key] = config[section][key].format(ENVDIR=env_dir, PIPEHOME=pipe_home)

    return config


def run_info(fname: str, reference: str, conda_env: str = "bp") -> None:
    config = read_config(reference, conda_env)
    pathlib.Path(os.path.dirname(fname)).mkdir(parents=True, exist_ok=True)
    with open(fname, "w") as run_file:
        run_file.write(
            "#PATH\nPIPE_HOME={}\nENV_DIR={}\n".format(
                config["PATH"]["pipe_home"], config["PATH"]["env_dir"]
            )
        )
        for section in ["TOOLS", "RESOURCES"]:
            run_file.write(f"\n#{section}\n")
            for key in config[section]:
                run_file.write(f"{key.upper()}={config[section][key]}\n")


def run_info_append(fname: str, line: str) -> None:
    with open(fname, "a") as run_file:
        run_file.write(line + "\n")


def log_dir(sample: str) -> str:
    log_dir = sample + "/logs"
    pathlib.Path(log_dir).mkdir(parents=True, exist_ok=True)
    return log_dir


def save_hold_jid(fname: str, jid: str) -> None:
    os.makedirs(os.path.dirname(fname), exist_ok=True)
    with open(fname, "w") as f:
        print(jid, file=f)
