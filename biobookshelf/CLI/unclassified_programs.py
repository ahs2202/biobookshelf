from biobookshelf.main import *
from biobookshelf import PKG


# Command line programs
def Recursively_stop_subprocesses(int_pid, df_processes=None):
    if int_pid < 1e3:
        print(
            "ERROR: given PID {} has less than 4 digits. For security reasons, this program failed".format(
                int_pid
            )
        ), sys.exit(2)
    if df_processes is None:  # get list of currently running processes
        df_processes = OS_Currently_running_processes()
    if int_pid in df_processes.PID.values:
        # print( f'terminating {int_pid}' )
        os.system("kill {}".format(int_pid))
    df_subprocesses = PD_Select(df_processes, PPID=int_pid)
    for int_pid_subprocess in df_subprocesses.PID.values:
        Recursively_stop_subprocesses(int_pid_subprocess, df_processes=df_processes)


def Stop_a_job(pid=None):
    """
    # 2022/04/05
    Stop a process with a given PID and all subprocess belonging to the PID.

    'pid' : pid to recursively terminates
    """
    flag_entry_point = PKG.Detect_Entry_Point(
        "biobook"
    )  # detect whether an entry point was used
    if flag_entry_point:
        parser = argparse.ArgumentParser(
            description="Stop a process with a given PID and all subprocess belonging to the PID. This program has been developed by Hyunsu An (2021/06/03)."
        )

        parser.add_argument(
            "-i", "--pid", help="PID of a job to be terminated. (required)", type=int
        )
        args = parser.parse_args()

        # [input] parse arguments from parse_args
        int_pid = args.pid

    """ [parse arguments] """
    if int_pid is None:
        print("required arguments are not given, exiting")
        if flag_entry_point:
            sys.exit()
        else:
            return -1
    Recursively_stop_subprocesses(int_pid)  # stop a process


def Server_Status():
    """
    # 2022/04/05
    Print CPU and MEMORY usage of the FGL server for each program for each user.
    """
    flag_entry_point = PKG.Detect_Entry_Point(
        "biobook"
    )  # detect whether an entry point was used
    if flag_entry_point:
        parser = argparse.ArgumentParser(
            description="Print CPU and MEMORY usage of the FGL server for each program for each user. This program has been developed by Hyunsu An."
        )
        args = parser.parse_args()

    df = Parse_Printed_Table(
        "\n".join(os.popen("top -bn 1").read().split("\n")[6:])
    ).rename(
        columns={"%CPU": "CPU", "%MEM": "MEM"}
    )  # rename columns to make it compatible with pandas DataFrame
    df = df[
        (df.CPU > 0) | (df.MEM > 0)
    ]  # if remove processes using '0' CPU and '0' MEM.
    print(" ----- Usage Total ----- \n")
    print(df[["CPU", "MEM"]].sum())
    print("\n\n ----- Usage by each user ----- \n")
    print(
        df[["USER", "CPU", "MEM"]]
        .groupby("USER")
        .sum()
        .sort_values("MEM", ascending=False)
    )
    print("\n\n ----- Usage by each program of each user ----- \n")
    print(
        df[["USER", "COMMAND", "CPU", "MEM"]]
        .groupby(["USER", "COMMAND"])
        .sum()
        .sort_values("MEM", ascending=False)
    )
    print()


def Release_Memory_from_Orphaned_Idle_Processes(
    *l_query, flag_dry_run: bool = False, flag_match_at_least_one_query: bool = True
):
    """
    * l_query, # list of query strings
    flag_dry_run : bool = False, # if True, does not terminate the processes and simply return the list of processes that will be terminated.
    flag_match_at_least_one_query : bool = True, # if True, only single query match will cause the process to be terminated. if False, the CMD line of the process should be matched with all queries in 'l_query'
    # 2024-01-09 12:44:35
    """
    flag_entry_point = PKG.Detect_Entry_Point(
        "biobook"
    )  # detect whether an entry point was used
    if flag_entry_point:
        parser = argparse.ArgumentParser(
            description="Release Orphaned Idle Processes (Interactive Python and R Jupyter kernels). This program has been developed by Hyunsu An."
        )
        args = parser.parse_args()

    # set default queries
    if len(l_query) == 0:
        l_query = [
            "python -m ipykernel_launcher -f ",
            "R --slave -e IRkernel::main() --args ",
        ]

    # get username of the current user
    user_name = os.getlogin()

    df_process = OS_Currently_running_processes()
    df_process = PD_Select(df_process, UID=user_name, PPID=1)
    name_col_cmd = list({"CMD", "TIME CMD"}.intersection(df_process.columns.values))[
        0
    ]  # some process brakes bk.Parse_Printed_Table by printing showing output that is not matched with the table format. in that case, TIME and CMD columns are combined into a single column. # identify the column name containing the CMD line
    arr_cmd_lines = df_process[name_col_cmd].values

    if flag_match_at_least_one_query:
        mask = np.zeros(len(arr_cmd_lines), dtype=bool)
        for query in l_query:
            mask |= Search_list_of_strings_with_multiple_query(
                arr_cmd_lines, query, return_mask=True
            )
    else:
        mask = Search_list_of_strings_with_multiple_query(
            arr_cmd_lines, *l_query, return_mask=True
        )
    df_process = df_process.loc[
        mask
    ]  # search the processes matched with at least one of the search queries

    print(f"{len( df_process )} number of processes will be terminated.")

    if not flag_dry_run:  # release memory of the matched orphaned processes
        for int_pid in df_process.PID.values:
            os.system(f"kill {int_pid}")

    return df_process
