#!/usr/bin/env python3

import subprocess as sp
import time
import getpass

NAME_MAP = {
    "bo791269": "Sichen",
    "gs418756": "Georg/MC",
    "nz257662": "Niko",
    "rs429310": "Robin",
    "sk656163": "Simon",
    "yx505956": "Yu Xu"
}

COLOR_FORMAT = "\033[{}m"
COLOR_RESET = COLOR_FORMAT.format("0")
COLOR_RED = COLOR_FORMAT.format("31")
COLOR_GREEN = COLOR_FORMAT.format("32")
COLOR_BLUE = COLOR_FORMAT.format("34")
COLOR_BOLD = COLOR_FORMAT.format("1")

def print_at(row, column, text, *args, **kwargs):
    print(f"\x1b7\x1b[{row};{column}f{text}\x1b8", *args, **kwargs)

def clear_terminal():
    print("\x1b[2J\x1b7\x1b[1;1f", end="")


class User:
    def __init__(self, name, running=0, total=0):
        self.name = name
        self.running = running
        self.total = total

    def __repr__(self):
        return f"{self.name} ({self.running}/{self.total})"


def get_jobs(project=None):
    command = ["bjobs", "-u", "all", "-o", "user stat user_group queue proj"]
    if project is not None:
        command.extend(["-P", project])
    text = sp.run(command, stdout=sp.PIPE, universal_newlines=True).stdout
    lines = [line.split(" ") for line in text.splitlines()]
    header = lines[0]
    jobs = [
        {
            key: value
            for key, value in zip(header, line)
        }
        for line in lines[1:]
    ]
    return jobs


def group_jobs(jobs):
    users = {}
    total_user = User("Total")
    for job in jobs:
        user_name = job["USER"]
        user = users.get(user_name, None) or User(user_name)
        users[user_name] = user
        if job["STAT"] == "RUN":
            user.running += 1
            total_user.running += 1
        user.total += 1
        total_user.total += 1
    return users, total_user


def format_user(user, color=False):
    name = NAME_MAP.get(user.name, user.name)
    name_str = name.ljust(8)
    run_str = str(user.running).rjust(4)
    total_str = str(user.total).rjust(4)
    if color:
        if getpass.getuser() == user.name:
            name_str = COLOR_BOLD + name_str + COLOR_RESET
        if user.total > 0:
            if user.running == 0:
                run_str = COLOR_RED + run_str + COLOR_RESET
            elif user.running < user.total:
                run_str = COLOR_BLUE + run_str + COLOR_RESET
            elif user.running == user.total:
                run_str = COLOR_GREEN + run_str + COLOR_RESET
    return f"{name_str} {run_str}/{total_str}"


def update(arguments):
    jobs = get_jobs(project=arguments.project)
    users, total_user = group_jobs(jobs)
    lines = []
    for row, name in enumerate(sorted(users)):
        user = users[name]
        print_name = NAME_MAP.get(name, name)
        lines.append(format_user(user, color=arguments.color))
    lines.append(format_user(total_user, color=arguments.color))
    print("\n".join(lines), end="")


def watch(arguments):
    while True:
        clear_terminal()
        update(arguments)
        time.sleep(arguments.delta)

def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--project", "-P", choices=["jara0052", "thes0375"])
    parser.add_argument("--watch", "-w", action="store_true")
    parser.add_argument("--delta", "-n", type=int, default=10)
    parser.add_argument("--color", "-c", action="store_true")

    arguments = parser.parse_args()
    if arguments.watch:
        watch(arguments)
    else:
        update(arguments)

if __name__ == "__main__":
    main()
