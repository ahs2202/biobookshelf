from biobookshelf.main import *
from biobookshelf import MAP
import biobookshelf as bk

pd.options.mode.chained_assignment = None  # default='warn' # to disable worining
import scanpy

from typing import Union, List, Literal, Dict, Callable, Set, Iterable, Tuple
import os
import pandas as pd
import numpy as np
import multiprocessing as mp
import math
import logging
from copy import deepcopy
import pickle
import time
import glob
import gzip  # to handle gzip file
import shutil  # for copying file
import base64  # for converting binary to text data (web application)
import json  # to read and write JSON file
import matplotlib.pyplot as plt
import scipy.sparse
import io
import concurrent.futures  # for multiprocessing


def is_s3_url(url):
    """# 2022-12-02 18:23:18
    check whether the given url is s3uri (s3url)
    """
    # handle None value
    if url is None:
        return False
    return "s3://" == url[:5]


def is_http_url(url):
    """# 2022-12-02 18:23:18
    check whether the given url is HTTP URL
    """
    return "https://" == url[:8] or "http://" == url[:7]


def is_remote_url(url):
    """# 2022-12-02 18:31:45
    check whether a url is a remote resource
    """
    return is_s3_url(url) or is_http_url(url)


""" remote files over HTTP """


def http_response_code(url):
    """# 2022-08-05 22:27:27
    check http response code
    """
    import requests  # download from url

    status_code = None  # by default, 'status_code' is None
    try:
        r = requests.head(url)
        status_code = r.status_code  # record the status header
    except requests.ConnectionError:
        status_code = None
    return status_code


def http_download_file(url, path_file_local):
    """# 2022-08-05 22:14:30
    download file from the remote location to the local directory
    """
    import requests  # download from url

    with requests.get(url, stream=True) as r:
        with open(path_file_local, "wb") as f:
            shutil.copyfileobj(r.raw, f)


""" remote files over AWS S3 """


def s3_exists(s3url):
    """# 2022-12-02 18:15:49
    check whether a path/file exists in AWS S3
    """
    import s3fs

    fs = s3fs.S3FileSystem()
    return fs.exists(s3url)


def s3_download_file(s3url, path_file_local):
    """# 2022-12-02 18:15:44
    download file from the remote AWS S3 location to the local directory
    """
    import s3fs

    fs = s3fs.S3FileSystem()
    fs.download(s3url, path_file_local)


def s3_rm(s3url, recursive=False, **kwargs):
    """# 2022-12-03 23:48:26
    delete file (or an entire folder) from a AWS S3 location
    """
    import s3fs

    fs = s3fs.S3FileSystem()
    fs.rm(s3url, recursive=recursive, **kwargs)  # delete files


""" method and class for handling file system """


# functions for various file system access
def filesystem_operations(
    method: Literal["exists", "rm", "glob", "mkdir", "mv", "cp", "isdir"],
    path_src: str,
    path_dest: Union[str, None] = None,
    flag_recursive: bool = True,
    dict_kwargs_credentials_s3: dict = dict(),
    **kwargs,
):
    """# 2022-12-04 00:57:45
    perform a file system operation (either Amazon S3 or local file system)

    method : Literal[
        'exists', # check whether a file or folder exists, given through 'path_src' arguments
        'rm', # remove file or folder, given through 'path_src' arguments
        'glob', # retrieve path of files matching the glob pattern, given through 'path_src' arguments
        'mkdir', # create a directory, given through 'path_src' arguments
        'mv', # move file or folder , given through 'path_src' and 'path_dest' arguments
        'cp', # copy file or folder , given through 'path_src' and 'path_dest' arguments
        'isdir', # check whether the given input is a file or directory
    ]

    kwargs :
        exist_ok : for 'mkdir' operation

    dict_kwargs_credentials_s3 : dict = dict( ) # the credentials for the Amazon S3 file system as keyworded arguments

    """
    if is_s3_url(path_src) or is_s3_url(
        path_dest
    ):  # if at least one path is s3 locations
        # %% Amazon s3 file system %%
        # load the file system
        import s3fs

        fs = s3fs.S3FileSystem(**dict_kwargs_credentials_s3)
        if method == "exists":
            return fs.exists(path_src, **kwargs)
        elif method == "rm":
            return fs.rm(path_src, recursive=flag_recursive, **kwargs)  # delete files
        elif method == "glob":
            return list(
                "s3://" + e for e in fs.glob(path_src, **kwargs)
            )  # 's3://' prefix should be added
        elif method == "mkdir":
            # use default 'exist_ok' value
            if "exist_ok" not in kwargs:
                kwargs["exist_ok"] = True
            return fs.makedirs(path_src, **kwargs)
        elif method == "mv":
            if not fs.exists(
                path_dest, **kwargs
            ):  # avoid overwriting of the existing file
                return fs.mv(path_src, path_dest, recursive=flag_recursive, **kwargs)
            else:
                return "destionation file already exists, exiting"
        elif method == "cp":
            if is_s3_url(path_src) and is_s3_url(path_dest):  # copy from s3 to s3
                return fs.copy(path_src, path_dest, recursive=flag_recursive, **kwargs)
            elif is_s3_url(path_src):  # copy from s3 to local
                return fs.get(path_src, path_dest, recursive=flag_recursive, **kwargs)
            elif is_s3_url(path_dest):  # copy from local to s3
                return fs.put(path_src, path_dest, recursive=flag_recursive, **kwargs)
        elif method == "isdir":
            return fs.isdir(path_src)
    elif is_http_url(path_src):  # for http
        # %% HTTP server %%
        if method == "exists":
            return (
                http_response_code(path_src) == 200
            )  # check whether http file (not tested for directory) exists
        else:
            return "not implemented"
    else:
        # %% local file system %%
        if method == "exists":
            return os.path.exists(path_src)
        elif method == "rm":
            if flag_recursive and os.path.isdir(
                path_src
            ):  # when the recursive option is active
                shutil.rmtree(path_src)
            else:
                os.remove(path_src)
        elif method == "glob":
            return glob.glob(path_src)
        elif method == "mkdir":
            # use default 'exist_ok' value
            if "exist_ok" not in kwargs:
                kwargs["exist_ok"] = True
            os.makedirs(path_src, exist_ok=kwargs["exist_ok"])
        elif method == "mv":
            shutil.move(path_src, path_dest)
        elif method == "cp":
            if flag_recursive and os.path.isdir(
                path_src
            ):  # when the recursive option is active
                shutil.copytree(path_src, path_dest)
            else:
                shutil.copyfile(path_src, path_dest)
        elif method == "isdir":
            return os.path.isdir(path_src)


""" previosuly written for biobookshelf """


def CB_Parse_list_of_id_cell(l_id_cell, dropna=True):
    """# 2022-03-25 16:35:23
    parse a given list of id_cells into a dataframe using 'SC.CB_detect_cell_barcode_from_id_cell' function
    'dropna' : drop id_cells that does not contains cell barcodes
    """
    df = pd.DataFrame(
        list([e] + list(CB_detect_cell_barcode_from_id_cell(e)) for e in l_id_cell),
        columns=["id_cell", "CB", "id_sample_from_id_cell"],
    ).set_index("id_cell")
    return df


def CB_Build_dict_id_sample_to_set_cb(l_id_cell):
    """# 2022-03-28 22:24:30
    Build a set of cell barcode for each id_sample from the given list of id_cells
    """
    df = CB_Parse_list_of_id_cell(l_id_cell)
    dict_id_sample_to_set_cb = dict()
    for cb, id_sample in df.values:
        if id_sample not in dict_id_sample_to_set_cb:
            dict_id_sample_to_set_cb[id_sample] = set()
        dict_id_sample_to_set_cb[id_sample].add(cb)
    return dict_id_sample_to_set_cb


def CB_Match_Batch(
    l_id_cell_1, l_id_cell_2, flag_calculate_proportion_using_l_id_cell_2=True
):
    """# 2022-03-28 23:10:43
    Find matching batches between two given lists of id_cells by finding the batches sharing the largest number of cell barcodes

    'l_id_cell_1' : first list of id_cells (e.g. unannotated barcodes)
    'l_id_cell_2' : second list of id_cells (e.g. annotated barcodes)
    'flag_calculate_proportion_using_l_id_cell_2' : if True, finding matching batches using the shared proportion calculated using the cell barcodes from 'l_id_cell_2'. if False, proportion of the matching barcodes will be calculated using the cell barcodes from 'l_id_cell_1'

    return:
    df_id_cell_matched, df_sample_matched
    """
    # retrieve set of cell barcodes
    df_id_cell_1 = CB_Parse_list_of_id_cell(l_id_cell_1)
    df_id_cell_2 = CB_Parse_list_of_id_cell(l_id_cell_2)
    dict_id_sample_to_set_cb_1 = CB_Build_dict_id_sample_to_set_cb(l_id_cell_1)
    dict_id_sample_to_set_cb_2 = CB_Build_dict_id_sample_to_set_cb(l_id_cell_2)

    # Find matching id_samples of the two given list of id_cells
    # calculate the proportion of matching cell barcodes between each pair of samples from the two given list of id_cells
    l_l = []
    for id_sample_1 in dict_id_sample_to_set_cb_1:
        for id_sample_2 in dict_id_sample_to_set_cb_2:
            set_cb_1 = dict_id_sample_to_set_cb_1[id_sample_1]
            set_cb_2 = dict_id_sample_to_set_cb_2[id_sample_2]
            float_prop_matching_cb = len(set_cb_1.intersection(set_cb_2)) / (
                len(set_cb_2)
                if flag_calculate_proportion_using_l_id_cell_2
                else len(set_cb_1)
            )
            l_l.append([id_sample_1, id_sample_2, float_prop_matching_cb])
    df = pd.DataFrame(
        l_l, columns=["id_sample_1", "id_sample_2", "float_prop_matching_cb"]
    )  # search result
    df_sample_matched = (
        df.sort_values("float_prop_matching_cb", ascending=False)
        .drop_duplicates("id_sample_2", keep="first")
        .drop_duplicates("id_sample_1", keep="first")
    )  # retrieve the best matches between samples so that a unique mapping exists for every sample

    # Find matching id_cells of given two list of id_cells
    df_id_cell_1.reset_index(inplace=True, drop=False)
    df_id_cell_2.reset_index(inplace=True, drop=False)
    df_id_cell_1.rename(
        columns={"id_sample_from_id_cell": "id_sample_from_id_cell_1"}, inplace=True
    )
    df_id_cell_2.rename(
        columns={"id_sample_from_id_cell": "id_sample_from_id_cell_2"}, inplace=True
    )
    df_id_cell_1["id_sample_from_id_cell_2"] = (
        df_id_cell_1.id_sample_from_id_cell_1.apply(
            MAP.Map(
                df_sample_matched.set_index("id_sample_1").id_sample_2.to_dict()
            ).a2b
        )
    )
    df_id_cell_1.dropna(
        subset=["id_sample_from_id_cell_2"], inplace=True
    )  # ignore cells without matching id_sample from the '2' batch
    df_id_cell_1.set_index(["CB", "id_sample_from_id_cell_2"], inplace=True)
    df_id_cell_matched = df_id_cell_1.join(
        df_id_cell_2[~pd.isnull(df_id_cell_2.id_sample_from_id_cell_2)].set_index(
            ["CB", "id_sample_from_id_cell_2"]
        ),
        lsuffix="_1",
        rsuffix="_2",
    )  # match id_cells from two given list of id_cells
    df_id_cell_matched.reset_index(drop=False, inplace=True)
    df_id_cell_matched = df_id_cell_matched[
        [
            "id_cell_1",
            "id_cell_2",
            "CB",
            "id_sample_from_id_cell_1",
            "id_sample_from_id_cell_2",
        ]
    ]  # reorder columns

    return df_id_cell_matched, df_sample_matched


def SC_Add_metadata(adata, df_metadata, suffix: str = "", inplace: bool = True):
    """# 2023-09-08 14:20:58
    Add metadata to the given AnnData

    adata
    df_metadata
    suffix : str = '', the suffix to the columns that will be added to AnnData.obs
    inplace : bool = True #
    """
    df_matched, _ = CB_Match_Batch(adata.obs.index.values, df_metadata.index.values)

    df_metadata = df_metadata.copy()  # copy the metadata dataframe
    df_matched.dropna(subset=["id_cell_2", "id_cell_1"], inplace=True)
    df_metadata["id_cell_1"] = df_matched.set_index("id_cell_2")[
        "id_cell_1"
    ]  # map id_cell
    df_metadata.dropna(subset=["id_cell_1"], inplace=True)  # drop invalid cells
    df_metadata.set_index("id_cell_1", inplace=True)
    if not inplace:
        adata = adata.copy()  # copy the anndata
    adata.obs = adata.obs.join(df_metadata, rsuffix=suffix)  # add the metadata

    return adata  # return the anndata


def SCANPY_Detect_cell_barcode_from_cell_id(adata):
    """# 2022-03-24 20:35:22
    Detect cell barcodes from id_cell (index of adata.obs), and add new two columns to the adata.obs [ 'CB', 'id_sample_from_id_cell' ]
    """
    adata.obs = adata.obs.join(
        pd.DataFrame(
            list(
                [e] + list(CB_detect_cell_barcode_from_id_cell(e))
                for e in adata.obs.index.values
            ),
            columns=["id_cell", "CB", "id_sample_from_id_cell"],
        ).set_index("id_cell")
    )


def SCANPY_Retrieve_Markers_as_DataFrame(adata):
    """# 2022-02-15 14:40:02
    receive scanpy anndata and return a dataframe contianing marker genes

    --- return ---
    df_marker : a dataframe contianing marker genes
    """
    l_df = []
    for index_clus, name_clus in enumerate(
        adata.uns["rank_genes_groups"]["names"].dtype.names
    ):
        df = pd.DataFrame(
            dict(
                (name_col, adata.uns["rank_genes_groups"][name_col][name_clus])
                for name_col in [
                    "logfoldchanges",
                    "names",
                    "pvals",
                    "pvals_adj",
                    "scores",
                ]
            )
        )
        df["name_clus"] = name_clus
        df["index_clus"] = index_clus
        l_df.append(df)
    df_marker = pd.concat(l_df)
    return df_marker


def CB_detect_cell_barcode_from_id_cell(
    id_cell, int_min_number_atgc_in_cell_barcode=16
):
    """# 2023-04-02 14:10:46
    retrieve cell_barcode from id_cell
    'int_min_number_atgc_in_cell_barcode' : number of ATGC characters in the cell barcode
    """
    int_count_atgc = 0
    int_start_appearance_of_atgc = None
    set_atgc = set("ATGC")

    def __retrieve_cell_barcode_and_id_channel_from_id_cell__(
        id_cell, int_start_appearance_of_atgc, int_count_atgc
    ):
        """__retrieve_cell_barcode_and_id_channel_from_id_cell__"""
        int_cb_start = int_start_appearance_of_atgc
        int_cb_end = int_start_appearance_of_atgc + int_count_atgc
        return [
            id_cell[int_cb_start:int_cb_end],
            id_cell[:int_cb_start] + "|" + id_cell[int_cb_end:],
        ]  # return cell_barcode, id_channel

    for index_c, c in enumerate(
        id_cell.upper()
    ):  # case-insensitive detection of cell-barcodes
        if c in set_atgc:
            if int_start_appearance_of_atgc is None:
                int_start_appearance_of_atgc = index_c
            int_count_atgc += 1
        else:
            """identify cell barcode and return the cell barcode"""
            if int_start_appearance_of_atgc is not None:
                if int_count_atgc >= int_min_number_atgc_in_cell_barcode:
                    return __retrieve_cell_barcode_and_id_channel_from_id_cell__(
                        id_cell, int_start_appearance_of_atgc, int_count_atgc
                    )
            # initialize the next search
            int_count_atgc = 0
            int_start_appearance_of_atgc = None
    """ identify cell barcode and return the cell barcode """
    if int_start_appearance_of_atgc is not None:
        if int_count_atgc >= int_min_number_atgc_in_cell_barcode:
            return __retrieve_cell_barcode_and_id_channel_from_id_cell__(
                id_cell, int_start_appearance_of_atgc, int_count_atgc
            )
    """ return None when cell_barcode was not found """
    return [None, None]


def Read_10X(path_folder_mtx_10x, verbose=False):
    """# 2021-11-24 13:00:13
    read 10x count matrix
    'path_folder_mtx_10x' : a folder containing files for 10x count matrix
    return df_mtx, df_feature
    """
    # handle inputs
    if path_folder_mtx_10x[-1] != "/":
        path_folder_mtx_10x += "/"

    # define input file directories
    path_file_bc = f"{path_folder_mtx_10x}barcodes.tsv.gz"
    path_file_feature = f"{path_folder_mtx_10x}features.tsv.gz"
    path_file_mtx = f"{path_folder_mtx_10x}matrix.mtx.gz"

    # check whether all required files are present
    if sum(
        list(
            not filesystem_operations("exists", path_folder)
            for path_folder in [path_file_bc, path_file_feature, path_file_mtx]
        )
    ):
        if verbose:
            logger.info(f"required file(s) is not present at {path_folder_mtx_10x}")

    # read mtx file as a tabular format
    df_mtx = pd.read_csv(path_file_mtx, sep=" ", comment="%")
    df_mtx.columns = ["id_row", "id_column", "read_count"]

    # read barcode and feature information
    df_bc = pd.read_csv(path_file_bc, sep="\t", header=None)
    df_bc.columns = ["barcode"]
    df_feature = pd.read_csv(path_file_feature, sep="\t", header=None)
    df_feature.columns = ["id_feature", "feature", "feature_type"]

    # mapping using 1 based coordinates (0->1 based coordinate )
    df_mtx["barcode"] = df_mtx.id_column.apply(
        bk.MAP.Map(
            bk.DICTIONARY_Build_from_arr(df_bc.barcode.values, index_start=1)
        ).a2b
    )  # mapping using 1 based coordinates (0->1 based coordinate )
    df_mtx["id_feature"] = df_mtx.id_row.apply(
        bk.Map(
            bk.DICTIONARY_Build_from_arr(df_feature.id_feature.values, index_start=1)
        ).a2b
    )
    df_mtx.drop(
        columns=["id_row", "id_column"], inplace=True
    )  # drop unnecessary columns

    return df_mtx, df_feature


def Write_10X(df_mtx, df_feature, path_folder_output_mtx_10x):
    """# 2021-11-24 12:57:30
    'df_feature' should contains the following column names : [ 'id_feature', 'feature', 'feature_type' ]
    'df_mtx' should contains the following column names : [ 'id_feature', 'barcode', 'read_count' ]
    'path_folder_output_mtx_10x' : an output folder directory where the mtx_10x files will be written

    """
    import scipy.io

    df_mtx = deepcopy(df_mtx)  # create a copy of df_mtx before modification

    # create an output folder
    filesystem_operations("mkdir", path_folder_output_mtx_10x, exist_ok=True)

    """ save barcode file """
    # retrieve list of barcodes
    arr_barcode = bk.LIST_COUNT(df_mtx.barcode, duplicate_filter=None).index.values
    pd.DataFrame(arr_barcode).to_csv(
        f"{path_folder_output_mtx_10x}barcodes.tsv.gz",
        sep="\t",
        index=False,
        header=False,
    )

    """ save feature file """
    # compose a feature dataframe
    df_feature[["id_feature", "feature", "feature_type"]].to_csv(
        f"{path_folder_output_mtx_10x}features.tsv.gz",
        sep="\t",
        index=False,
        header=False,
    )  # save as a file
    # retrieve list of features
    arr_id_feature = df_feature.id_feature.values

    """ save matrix file """
    # convert feature and barcode to integer indices
    df_mtx.id_feature = df_mtx.id_feature.apply(
        bk.Map(
            bk.DICTIONARY_Build_from_arr(arr_id_feature, order_index_entry=False)
        ).a2b
    )  # 0-based coordinates
    df_mtx.barcode = df_mtx.barcode.apply(
        bk.Map(bk.DICTIONARY_Build_from_arr(arr_barcode, order_index_entry=False)).a2b
    )  # 0-based coordinates
    # save count matrix as a gzipped matrix market format
    row, col, data = df_mtx[["id_feature", "barcode", "read_count"]].values.T
    sm = scipy.sparse.coo_matrix(
        (data, (row, col)), shape=(len(arr_id_feature), len(arr_barcode))
    )
    scipy.io.mmwrite(f"{path_folder_output_mtx_10x}matrix", sm)
    # remove previous output file to overwrite the file
    path_file_mtx_output = f"{path_folder_output_mtx_10x}matrix.mtx.gz"
    if filesystem_operations("exists", path_file_mtx_output):
        filesystem_operations("rm", path_file_mtx_output)
    bk.OS_Run(["gzip", f"{path_folder_output_mtx_10x}matrix.mtx"])  # gzip the mtx file


def AnnData_Convert_to_10X_MTX(
    adata,
    path_folder_mtx_output,
    dict_var_rename: dict = {"feature_types": "feature_type", "gene_ids": "id_feature"},
    dtype_value=np.int64,
):
    """# 2022-12-14 02:14:31
    write AnnData count matrix as a 10X matrix object

    'dict_var_rename' : a dictionary for renaming columns of adata.var columns
    """
    import scipy.io

    # compose df_var
    df_feature = adata.var
    df_feature.rename(columns=dict_var_rename, inplace=True)

    # create an output folder
    filesystem_operations("mkdir", path_folder_mtx_output, exist_ok=True)

    """ save barcode file """
    # retrieve list of barcodes
    arr_barcode = adata.obs.index.values
    pd.DataFrame(arr_barcode).to_csv(
        f"{path_folder_mtx_output}barcodes.tsv.gz", sep="\t", index=False, header=False
    )

    """ save feature file """
    # compose a feature dataframe
    df_feature[["id_feature", "feature", "feature_type"]].to_csv(
        f"{path_folder_mtx_output}features.tsv.gz", sep="\t", index=False, header=False
    )  # save as a file
    # retrieve list of features
    arr_id_feature = df_feature.id_feature.values

    """ save matrix file """
    # save count matrix as a gzipped matrix market format
    arr_int_barcode, arr_int_id_feature, arr_read_count = scipy.sparse.find(adata.X)
    # convert dtype of the values
    if dtype_value is not None:
        arr_read_count = arr_read_count.astype(dtype_value)
    # compose a sparse matrix
    sm = scipy.sparse.coo_matrix(
        (arr_read_count, (arr_int_id_feature, arr_int_barcode)),
        shape=(len(arr_id_feature), len(arr_barcode)),
    )
    scipy.io.mmwrite(f"{path_folder_mtx_output}matrix", sm)
    # remove previous output file to overwrite the file
    path_file_mtx_output = f"{path_folder_mtx_output}matrix.mtx.gz"
    if filesystem_operations("exists", path_file_mtx_output):
        filesystem_operations("rm", path_file_mtx_output)
    bk.OS_Run(["gzip", f"{path_folder_mtx_output}matrix.mtx"])  # gzip the mtx file


def __function_for_adjusting_thresholds_for_filtering_empty_droplets__(
    path_folder_mtx_10x_output, min_counts, min_features, min_cells
):
    """# 2022-02-23 14:26:07
    This function is intended for the use in 'MTX_10X_Filter' function for filtering cells from the 10X dataset (before chromium X, 10,000 cells per channel)

    Assuming a typical number of droplets in a experiment is 100,000, adjust 'min_counts' to reduce the number of filtered cells below 'int_max_num_cells'
    """
    s_count = (
        pd.read_csv(
            f"{path_folder_mtx_10x_output}dict_id_column_to_count.before_filtering.tsv.gz",
            sep="\t",
            header=None,
            index_col=0,
        )[1]
        .sort_values(ascending=False)
        .iloc[:100000]
    )

    int_max_num_cells = 20000  # maximum number of allowed cells
    min_counts_maximum = 2000

    def function_for_increasing_min_counts(min_counts):
        return min_counts * 2

    while True:
        """increase threshold if the number of filtered cells is larger than 'int_max_num_cells'"""
        if (
            len(s_count[s_count > min_counts]) > int_max_num_cells
            and min_counts < min_counts_maximum
        ):
            min_counts = function_for_increasing_min_counts(min_counts)
        else:
            break
    return min_counts, min_features, min_cells


def _get_matrix_market_data_type_from_line(line):
    """# 2023-09-10 16:06:24
    analyze the line to get matrix market data type in string format ('real' or 'integer')
    """
    try:
        return line.split("%%MatrixMarket")[1].strip().split()[2]
    except:
        return  # if an error occurs, return None


def _MTX_Detect_data_type(path_file_mtx):
    """# 2023-09-10 16:01:34
    detect data type of the input matrix file
    """
    with (
        gzip.open(path_file_mtx, "rt")
        if ".gz" == path_file_mtx[-3:]
        else open(path_file_mtx, "r")
    ) as file:
        str_datatype = _get_matrix_market_data_type_from_line(
            file.readline()
        )  # matrix market header line is present in the first line
    return str_datatype  # return the datatype


def MTX_10X_Split(
    path_folder_mtx_10x_output,
    int_max_num_entries_for_chunk=10000000,
    flag_split_mtx=True,
    flag_split_mtx_again=False,
):
    """# 2022-04-28 01:16:35
    split input mtx file into multiple files and write a flag file indicating the splitting has been completed.
    return the list of split mtx files

    'flag_split_mtx' : if 'flag_split_mtx' is True, split input mtx file into multiple files. if False, does not split the input matrix, and just return the list containing a single path pointing to the input matrix. This flag exists for the compatibility with single-thread operations
    'flag_split_mtx_again' : split the input matrix again even if it has beem already split. It will remove previously split files.
    """
    # 'flag_split_mtx' : if False, does not split the input matrix, and just return the list containing a single path pointing to the input matrix
    if not flag_split_mtx:
        return [f"{path_folder_mtx_10x_output}matrix.mtx.gz"]

    """ if 'flag_split_mtx_again' flag is True, remove previously split files """
    path_file_flag = f"{path_folder_mtx_10x_output}matrix.mtx.gz.split.flag"
    if flag_split_mtx_again:
        filesystem_operations("rm", path_file_flag)  # remove the flag
        # remove previously split files
        for path_file in filesystem_operations(
            "glob", f"{path_folder_mtx_10x_output}matrix.mtx.gz.*.gz"
        ):
            filesystem_operations("rm", path_file)

    """ split input matrix file """
    if not filesystem_operations(
        "exists", path_file_flag
    ):  # check whether the flag exists
        index_mtx_10x = 0
        newfile = gzip.open(
            f"{path_folder_mtx_10x_output}matrix.mtx.gz.{index_mtx_10x}.gz", "wb"
        )
        l_path_file_mtx_10x = [
            f"{path_folder_mtx_10x_output}matrix.mtx.gz.{index_mtx_10x}.gz"
        ]
        int_num_entries_written_for_the_current_chunk = 0
        with gzip.open(f"{path_folder_mtx_10x_output}matrix.mtx.gz", "rb") as file:
            while True:
                line = file.readline()  # binary string
                if len(line) == 0:
                    newfile.close()  # close the output file
                    break
                """ write the line to the current chunk and update the number of entries written for the current chunk """
                newfile.write(line)
                int_num_entries_written_for_the_current_chunk += 1
                """ initialize the next chunk if a sufficient number of entries were written """
                if (
                    int_num_entries_written_for_the_current_chunk
                    >= int_max_num_entries_for_chunk
                ):
                    newfile.close()  # close the output file
                    # initialize the next chunk
                    index_mtx_10x += 1
                    int_num_entries_written_for_the_current_chunk = 0
                    newfile = gzip.open(
                        f"{path_folder_mtx_10x_output}matrix.mtx.gz.{index_mtx_10x}.gz",
                        "wb",
                    )
                    l_path_file_mtx_10x.append(
                        f"{path_folder_mtx_10x_output}matrix.mtx.gz.{index_mtx_10x}.gz"
                    )
        with open(path_file_flag, "w") as file:
            file.write("completed")
    else:
        """retrieve the list of split mtx files"""
        df = bk.GLOB_Retrive_Strings_in_Wildcards(
            f"{path_folder_mtx_10x_output}matrix.mtx.gz.*.gz"
        )
        df.wildcard_0 = df.wildcard_0.astype(int)
        df.sort_values("wildcard_0", ascending=True, inplace=True)
        l_path_file_mtx_10x = df.path.values
    return l_path_file_mtx_10x


dict_id_feature_to_index_feature = dict()


def __MTX_10X_Combine__renumber_feature_mtx_10x__(
    path_file_input, path_folder_mtx_10x_output
):
    """# deprecated
    internal function for MTX_10X_Combine
    # 2022-02-22 00:38:33
    """
    #     dict_id_feature_to_index_feature = bk.PICKLE_Read( f'{path_folder_mtx_10x_output}dict_id_feature_to_index_feature.pickle' ) # retrieve id_feature to index_feature mapping
    for (
        path_folder_mtx_10x,
        int_total_n_barcodes_of_previously_written_matrices,
        index_mtx_10x,
    ) in pd.read_csv(path_file_input, sep="\t").values:
        # directly write matrix.mtx.gz file without header
        with gzip.open(
            f"{path_folder_mtx_10x_output}matrix.mtx.gz.{index_mtx_10x}.gz", "wb"
        ) as newfile:
            arr_id_feature = pd.read_csv(
                f"{path_folder_mtx_10x}features.tsv.gz", sep="\t", header=None
            ).values[
                :, 0
            ]  # retrieve a list of id_feature for the current dataset
            with gzip.open(
                f"{path_folder_mtx_10x}matrix.mtx.gz", "rb"
            ) as file:  # retrieve a list of features
                line = file.readline().decode()  # read the first line
                # if the first line of the file contains a comment line, read all comment lines and a description line following the comments.
                if len(line) > 0 and line[0] == "%":
                    # read comment and the description line
                    while True:
                        if line[0] != "%":
                            break
                        line = file.readline().decode()  # read the next line
                    line = (
                        file.readline().decode()
                    )  # discard the description line and read the next line
                # process entries
                while True:
                    if len(line) == 0:
                        break
                    index_row, index_col, int_value = tuple(
                        map(int, line.strip().split())
                    )  # parse each entry of the current matrix
                    newfile.write(
                        (
                            " ".join(
                                tuple(
                                    map(
                                        str,
                                        [
                                            dict_id_feature_to_index_feature[
                                                arr_id_feature[index_row - 1]
                                            ],
                                            index_col
                                            + int_total_n_barcodes_of_previously_written_matrices,
                                            int_value,
                                        ],
                                    )
                                )
                            )
                            + "\n"
                        ).encode()
                    )  # translate indices of the current matrix to that of the combined matrix
                    line = file.readline().decode()  # read the next line


def Read_SPLiT_Seq(path_folder_mtx):
    """# 2022-04-22 07:10:50
    Read SPLiT-Seq pipeline output
    return:
    df_feature, df_mtx
    """
    path_file_bc = f"{path_folder_mtx}cell_metadata.csv"
    path_file_feature = f"{path_folder_mtx}genes.csv"
    path_file_mtx = f"{path_folder_mtx}DGE.mtx"

    # read mtx file as a tabular format
    df_mtx = pd.read_csv(path_file_mtx, sep=" ", comment="%")
    df_mtx.columns = ["id_row", "id_column", "read_count"]

    # read barcode and feature information
    df_bc = pd.read_csv(path_file_bc)[["cell_barcode"]]
    df_bc.columns = ["barcode"]
    df_feature = pd.read_csv(path_file_feature, index_col=0)
    df_feature.columns = ["id_feature", "feature", "genome"]

    # mapping using 1 based coordinates (0->1 based coordinate )
    df_mtx["barcode"] = df_mtx.id_row.apply(
        bk.Map(bk.DICTIONARY_Build_from_arr(df_bc.barcode.values, index_start=1)).a2b
    )  # mapping using 1 based coordinates (0->1 based coordinate )
    df_mtx["id_feature"] = df_mtx.id_column.apply(
        bk.Map(
            bk.DICTIONARY_Build_from_arr(df_feature.id_feature.values, index_start=1)
        ).a2b
    )
    df_mtx.drop(
        columns=["id_row", "id_column"], inplace=True
    )  # drop unnecessary columns
    return df_feature, df_mtx


def MTX_10X_Barcode_add_prefix_or_suffix(
    path_file_barcodes_input,
    path_file_barcodes_output=None,
    barcode_prefix="",
    barcode_suffix="",
):
    """# 2022-05-13 17:54:13
    Add prefix or suffix to the 'barcode' of a given 'barcodes.tsv.gz' file
    'path_file_barcodes_output' : default: None. by default, the input 'path_file_barcodes_input' file will be overwritten with the modified barcodes
    """
    flag_replace_input_file = (
        path_file_barcodes_output is None
    )  # retrieve a flag indicating the replacement of original input file with modified input file
    if flag_replace_input_file:
        path_file_barcodes_output = (
            f"{path_file_barcodes_input}.writing.tsv.gz"  # set a temporary output file
        )
    newfile = gzip.open(path_file_barcodes_output, "wb")  # open an output file
    with gzip.open(path_file_barcodes_input, "rb") as file:
        while True:
            line = file.readline()
            if len(line) == 0:
                break
            barcode = line.decode().strip()  # parse a barcode
            barcode_new = (
                barcode_prefix + barcode + barcode_suffix
            )  # compose new barcode
            newfile.write((barcode_new + "\n").encode())  # write a new barcode
    newfile.close()  # close the output file
    # if the output file path was not given, replace the original file with modified file
    if flag_replace_input_file:
        filesystem_operations("rm", path_file_barcodes_input)
        filesystem_operations("mv", path_file_barcodes_output, path_file_barcodes_input)


def MTX_10X_Feature_add_prefix_or_suffix(
    path_file_features_input,
    path_file_features_output=None,
    id_feature_prefix="",
    id_feature_suffix="",
    name_feature_prefix="",
    name_feature_suffix="",
):
    """# 2022-05-13 17:54:17
    Add prefix or suffix to the id_feature and name_feature of a given 'features.tsv.gz' file
    'path_file_features_output' : default: None. by default, the input 'path_file_features_input' file will be overwritten with the modified features
    """
    flag_replace_input_file = (
        path_file_features_output is None
    )  # retrieve a flag indicating the replacement of original input file with modified input file
    if flag_replace_input_file:
        path_file_features_output = (
            f"{path_file_features_input}.writing.tsv.gz"  # set a temporary output file
        )
    newfile = gzip.open(path_file_features_output, "wb")  # open an output file
    with gzip.open(path_file_features_input, "rb") as file:
        while True:
            line = file.readline()
            if len(line) == 0:
                break
            id_feature, name_feature, type_feature = (
                line.decode().strip().split("\t")
            )  # parse a feature
            id_feature_new = (
                id_feature_prefix + id_feature + id_feature_suffix
            )  # compose new id_feature
            name_feature_new = (
                name_feature_prefix + name_feature + name_feature_suffix
            )  # compose new name_feature
            newfile.write(
                (
                    "\t".join([id_feature_new, name_feature_new, type_feature]) + "\n"
                ).encode()
            )  # write a new feature
    newfile.close()  # close the output file
    # if the output file path was not given, replace the original file with modified file
    if flag_replace_input_file:
        filesystem_operations("rm", path_file_features_input)
        filesystem_operations("mv", path_file_features_output, path_file_features_input)


def __MTX_10X_Combine__renumber_barcode_or_feature_index_mtx_10x__(
    path_file_input: str,
    path_folder_mtx_10x_output: str,
    flag_renumber_feature_index: bool,
    str_datatype: str,
):
    """
    internal function for MTX_10X_Combine
    # 2023-09-10 17:15:05

    flag_renumber_feature_index : bool : if True, assumes barcodes are not shared between matrices and renumber features only. If False, assumes features are not shared between matrices and renumber barcodes only.
    str_data_type : str # matrix market datatype in string format
    """
    global dict_id_entry_to_index_entry
    flag_matrix_contain_float_values = (
        str_data_type == "real"
    )  # retrieve a flag indicating the matrix is containing float values
    for (
        path_folder_mtx_10x,
        int_total_n_entries_of_previously_written_matrices,
        index_mtx_10x,
    ) in pd.read_csv(path_file_input, sep="\t").values:
        # directly write matrix.mtx.gz file without header
        with gzip.open(
            f"{path_folder_mtx_10x_output}matrix.mtx.gz.{index_mtx_10x}.gz", "wb"
        ) as newfile:
            arr_id_entry = pd.read_csv(
                f"{path_folder_mtx_10x}{'features' if flag_renumber_feature_index else 'barcodes'}.tsv.gz",
                sep="\t",
                header=None,
            ).values[
                :, 0
            ]  # retrieve a list of id_feature for the current dataset
            with gzip.open(
                f"{path_folder_mtx_10x}matrix.mtx.gz", "rb"
            ) as file:  # retrieve a list of features
                line = file.readline().decode()  # read the first line
                # if the first line of the file contains a comment line, read all comment lines and a description line following the comments.
                if len(line) > 0 and line[0] == "%":
                    # read comment and the description line
                    while True:
                        if line[0] != "%":
                            break
                        line = file.readline().decode()  # read the next line
                    line = (
                        file.readline().decode()
                    )  # discard the description line and read the next line
                # process entries
                while True:
                    if len(line) == 0:
                        break
                    """ parse a record """
                    l_str = (
                        line.strip().split()
                    )  # parse a record of a matrix-market format file
                    if flag_matrix_contain_float_values:  # parse float data type
                        id_row, id_column, value = (
                            int(l_str[0]),
                            int(l_str[1]),
                            float(l_str[2]),
                        )
                    else:  # parse integer data type
                        id_row, id_column, value = tuple(int(float(e)) for e in l_str)

                    newfile.write(
                        (
                            " ".join(
                                tuple(
                                    map(
                                        str,
                                        (
                                            [
                                                dict_id_entry_to_index_entry[
                                                    arr_id_entry[id_row - 1]
                                                ],
                                                id_column
                                                + int_total_n_entries_of_previously_written_matrices,
                                            ]
                                            if flag_renumber_feature_index
                                            else [
                                                id_row
                                                + int_total_n_entries_of_previously_written_matrices,
                                                dict_id_entry_to_index_entry[
                                                    arr_id_entry[id_column - 1]
                                                ],
                                            ]
                                        )
                                        + [value],
                                    )
                                )
                            )
                            + "\n"
                        ).encode()
                    )  # translate indices of the current matrix to that of the combined matrix
                    line = file.readline().decode()  # read the next line


def MTX_10X_Combine(
    path_folder_mtx_10x_output,
    *l_path_folder_mtx_10x_input,
    int_num_threads=15,
    flag_split_mtx=True,
    flag_split_mtx_again=False,
    int_max_num_entries_for_chunk=10000000,
    flag_low_memory_mode_because_there_is_no_shared_cell_between_mtxs=None,
    flag_low_memory_mode_because_there_is_no_shared_feature_between_mtxs=None,
    verbose=False,
):
    """
    # 2023-09-10 17:14:57
    Combine 10X count matrix files from the given list of folders and write combined output files to the given output folder 'path_folder_mtx_10x_output'
    If there are no shared cells between matrix files, a low-memory mode will be used. The output files will be simply combined since no count summing operation is needed. Only feature matrix will be loaded and updated in the memory.
    'id_feature' should be unique across all features. if id_feature is not unique, features with duplicated id_features will lead to combining of the features into a single feature (with combined counts/values).

    'int_num_threads' : number of threads to use when combining datasets. multiple threads will be utilized only when datasets does not share cells and thus can be safely concatanated.
    'flag_split_mtx' : split the resulting mtx file so that the contents in the output mtx file can be processed in parallel without ungzipping the mtx.gz file and spliting the file.
    'flag_low_memory_mode_because_there_is_no_shared_cell_between_mtxs' : a flag for entering low-memory mode when there is no shared cells between given input matrices. By default (when None is given), matrices will be examined and the flag will be set automatically by the program. To reduce running time and memory, this flag can be manually set by users. Explicitly setting this flag will dramatically reduce the memory consumption.
    'flag_low_memory_mode_because_there_is_no_shared_feature_between_mtxs' : a flag for entering low-memory mode when there is no shared features between given input matrices. By default (when None is given), matrices will be examined and the flag will be set automatically by the program. To reduce running time and memory, this flag can be manually set by users. Explicitly setting this flag will dramatically reduce the memory consumption.
    """

    # create an output folder
    filesystem_operations("mkdir", path_folder_mtx_10x_output, exist_ok=True)

    str_datatype = (
        "real"
        if sum(
            _MTX_Detect_data_type(f"{path_folder}matrix.mtx.gz") == "real"
            for path_folder in l_path_folder_mtx_10x_input
        )
        else "integer"
    )  # retrieve datatype of the operation. if at least one matrix contains real data type, treat all matrices as matrix containing real data

    if (
        not flag_low_memory_mode_because_there_is_no_shared_feature_between_mtxs
        and flag_low_memory_mode_because_there_is_no_shared_cell_between_mtxs is None
    ):
        """retrieve cell barcodes of all 10X matrices and check whether cell barcodes are not shared between matrices"""
        int_total_n_barcodes_of_previously_written_matrices = (
            0  # follow the number of barcodes that are previously written
        )
        l_int_total_n_barcodes_of_previously_written_matrices = (
            []
        )  # calculate the number of barcodes of the previous dataset in the combined mtx.
        set_barcode = set()  # update a set of unique barcodes
        for path_folder_mtx_10x in l_path_folder_mtx_10x_input:
            arr_barcode = (
                pd.read_csv(
                    f"{path_folder_mtx_10x}barcodes.tsv.gz", sep="\t", header=None
                )
                .squeeze("columns")
                .values
            )  # retrieve a list of features
            set_barcode.update(arr_barcode)  # update a set of barcodes
            l_int_total_n_barcodes_of_previously_written_matrices.append(
                int_total_n_barcodes_of_previously_written_matrices
            )
            int_total_n_barcodes_of_previously_written_matrices += len(
                arr_barcode
            )  # update the number of barcodes
        """ check whether there are shared cell barcodes between matrices and set a flag for entering a low-memory mode """
        flag_low_memory_mode_because_there_is_no_shared_cell_between_mtxs = (
            len(set_barcode) == int_total_n_barcodes_of_previously_written_matrices
        )  # update flag
    elif flag_low_memory_mode_because_there_is_no_shared_cell_between_mtxs:
        """retrieve cell barcodes of all 10X matrices and check whether cell barcodes are not shared between matrices"""
        int_total_n_barcodes_of_previously_written_matrices = (
            0  # follow the number of barcodes that are previously written
        )
        l_int_total_n_barcodes_of_previously_written_matrices = (
            []
        )  # calculate the number of barcodes of the previous dataset in the combined mtx.
        for path_folder_mtx_10x in l_path_folder_mtx_10x_input:
            l_int_total_n_barcodes_of_previously_written_matrices.append(
                int_total_n_barcodes_of_previously_written_matrices
            )
            int_total_n_barcodes_of_previously_written_matrices += len(
                pd.read_csv(
                    f"{path_folder_mtx_10x}barcodes.tsv.gz", sep="\t", header=None
                )
            )  # retrieve a list of barcodes and # update the number of barcodes

    if (
        not flag_low_memory_mode_because_there_is_no_shared_cell_between_mtxs
        and flag_low_memory_mode_because_there_is_no_shared_feature_between_mtxs is None
    ):
        """retrieve features of all 10X matrices and check whether features are not shared between matrices"""
        int_total_n_features_of_previously_written_matrices = (
            0  # follow the number of features that are previously written
        )
        l_int_total_n_features_of_previously_written_matrices = (
            []
        )  # calculate the number of features of the previous dataset in the combined mtx.
        set_feature = set()  # update a set of unique features
        for path_folder_mtx_10x in l_path_folder_mtx_10x_input:
            arr_feature = (
                pd.read_csv(
                    f"{path_folder_mtx_10x}features.tsv.gz",
                    sep="\t",
                    header=None,
                    usecols=[0],
                )
                .squeeze("columns")
                .values
            )  # retrieve a list of features
            set_feature.update(arr_feature)  # update a set of features
            l_int_total_n_features_of_previously_written_matrices.append(
                int_total_n_features_of_previously_written_matrices
            )
            int_total_n_features_of_previously_written_matrices += len(
                arr_feature
            )  # update the number of features
        """ check whether there are shared features between matrices and set a flag for entering a low-memory mode """
        flag_low_memory_mode_because_there_is_no_shared_feature_between_mtxs = (
            len(set_feature) == int_total_n_features_of_previously_written_matrices
        )  # update flag
    elif flag_low_memory_mode_because_there_is_no_shared_feature_between_mtxs:
        """retrieve features of all 10X matrices and check whether features are not shared between matrices"""
        int_total_n_features_of_previously_written_matrices = (
            0  # follow the number of features that are previously written
        )
        l_int_total_n_features_of_previously_written_matrices = (
            []
        )  # calculate the number of features of the previous dataset in the combined mtx.
        for path_folder_mtx_10x in l_path_folder_mtx_10x_input:
            l_int_total_n_features_of_previously_written_matrices.append(
                int_total_n_features_of_previously_written_matrices
            )
            int_total_n_features_of_previously_written_matrices += len(
                pd.read_csv(
                    f"{path_folder_mtx_10x}features.tsv.gz",
                    sep="\t",
                    header=None,
                    usecols=[0],
                )
            )  # retrieve a list of features and update the number of features

    flag_low_memory_mode = (
        flag_low_memory_mode_because_there_is_no_shared_cell_between_mtxs
        or flag_low_memory_mode_because_there_is_no_shared_feature_between_mtxs
    )  # retrieve flag for low-memory mode
    if flag_low_memory_mode:
        """low-memory mode"""
        flag_renumber_feature_index = flag_low_memory_mode_because_there_is_no_shared_cell_between_mtxs  # retrieve a flag for renumbering features
        if verbose:
            logger.info(
                f"entering low-memory mode, re-numbering {'features' if flag_renumber_feature_index else 'barcodes'} index because {'barcodes' if flag_renumber_feature_index else 'features'} are not shared across the matrices."
            )

        """ write a combined barcodes/features.tsv.gz - that are not shared between matrices """
        bk.OS_Run(
            ["cat"]
            + list(
                f"{path_folder_mtx_10x}{'barcodes' if flag_renumber_feature_index else 'features'}.tsv.gz"
                for path_folder_mtx_10x in l_path_folder_mtx_10x_input
            ),
            path_file_stdout=f"{path_folder_mtx_10x_output}{'barcodes' if flag_renumber_feature_index else 'features'}.tsv.gz",
            stdout_binary=True,
            return_output=False,
        )  # combine the files in order

        """ collect a set of unique entries and a list of entries for each 10X matrix """
        set_t_entry = set()  # update a set unique id_entry (either id_cell or id_entry)
        for path_folder_mtx_10x in l_path_folder_mtx_10x_input:
            set_t_entry.update(
                list(
                    map(
                        tuple,
                        pd.read_csv(
                            f"{path_folder_mtx_10x}{'features' if flag_renumber_feature_index else 'barcodes'}.tsv.gz",
                            sep="\t",
                            header=None,
                        ).values,
                    )
                )
            )  # update a set of feature tuples

        """ write a combined features/barcodes.tsv.gz - that are shared between matrices """
        l_t_entry = list(set_t_entry)  # convert set to list
        with gzip.open(
            f"{path_folder_mtx_10x_output}{'features' if flag_renumber_feature_index else 'barcodes'}.tsv.gz",
            "wb",
        ) as newfile:
            for t_entry in l_t_entry:
                newfile.write(("\t".join(t_entry) + "\n").encode())

        """ build a mapping of id_entry to index_entry, which will be consistent across datasets - for features/barcodes that are shared between matrices """
        global dict_id_entry_to_index_entry  # use global variable for multiprocessing
        dict_id_entry_to_index_entry = dict(
            (t_entry[0], index_entry + 1)
            for index_entry, t_entry in enumerate(l_t_entry)
        )  # 0>1 based index
        bk.PICKLE_Write(
            f"{path_folder_mtx_10x_output}dict_id_entry_to_index_entry.pickle",
            dict_id_entry_to_index_entry,
        )  # save id_feature to index_feature mapping as a pickle file

        """ collect the number of records for each 10X matrix """
        int_total_n_records = 0
        for path_folder_mtx_10x in l_path_folder_mtx_10x_input:
            with gzip.open(
                f"{path_folder_mtx_10x}matrix.mtx.gz", "rt"
            ) as file:  # retrieve a list of features
                line = file.readline()
                while line[0] == "%":
                    line = file.readline()
                int_total_n_records += int(
                    line.strip().split()[2]
                )  # update the total number of entries

        """ write a part of a combined matrix.mtx.gz for each dataset using multiple processes """
        # compose inputs for multiprocessing
        df_input = pd.DataFrame(
            {
                "path_folder_input_mtx_10x": l_path_folder_mtx_10x_input,
                "int_total_n_barcodes_of_previously_written_matrices": (
                    l_int_total_n_barcodes_of_previously_written_matrices
                    if flag_renumber_feature_index
                    else l_int_total_n_features_of_previously_written_matrices
                ),
                "index_mtx_10x": np.arange(
                    len(l_int_total_n_barcodes_of_previously_written_matrices)
                    if flag_renumber_feature_index
                    else len(l_int_total_n_features_of_previously_written_matrices)
                ),
            }
        )
        bk.Multiprocessing(
            df_input,
            __MTX_10X_Combine__renumber_barcode_or_feature_index_mtx_10x__,
            int_num_threads,
            global_arguments=[
                path_folder_mtx_10x_output,
                flag_renumber_feature_index,
                str_datatype,
            ],
        )
        #         filesystem_operations( 'rm', f'{path_folder_mtx_10x_output}dict_id_entry_to_index_entry.pickle' ) # remove pickle file

        """ combine parts and add the MTX file header to compose a combined mtx file """
        df_file = bk.GLOB_Retrive_Strings_in_Wildcards(
            f"{path_folder_mtx_10x_output}matrix.mtx.gz.*.gz"
        )
        df_file.wildcard_0 = df_file.wildcard_0.astype(int)
        df_file.sort_values("wildcard_0", inplace=True)

        # write header
        path_file_header = f"{path_folder_mtx_10x_output}matrix.mtx.header.txt.gz"
        with gzip.open(path_file_header, "wt") as newfile:
            newfile.write(
                f"%%MatrixMarket matrix coordinate {str_datatype} general\n%\n{len( l_t_entry ) if flag_renumber_feature_index else int_total_n_features_of_previously_written_matrices} {int_total_n_barcodes_of_previously_written_matrices if flag_renumber_feature_index else len( l_t_entry )} {int_total_n_records}\n"
            )
        bk.OS_Run(
            ["cat", path_file_header] + list(df_file.path.values),
            path_file_stdout=f"{path_folder_mtx_10x_output}matrix.mtx.gz",
            stdout_binary=True,
            return_output=False,
        )  # combine the output mtx files in the order

        if not flag_split_mtx:
            # delete temporary files if 'flag_split_mtx' is False
            for path_file in df_file.path.values:
                os.remove(path_file)

        # write a flag indicating that the current output directory contains split mtx files
        with open(f"{path_folder_mtx_10x_output}matrix.mtx.gz.split.flag", "w") as file:
            file.write("completed")
    else:
        """normal operation mode perfoming count merging operations"""
        l_df_mtx, l_df_feature = [], []
        for path_folder_mtx_10x in l_path_folder_mtx_10x_input:
            df_mtx, df_feature = Read_10X(path_folder_mtx_10x)
            l_df_mtx.append(df_mtx), l_df_feature.append(df_feature)

        # combine mtx
        df_mtx = pd.concat(l_df_mtx)
        df_mtx = df_mtx.groupby(["barcode", "id_feature"]).sum()
        df_mtx.reset_index(drop=False, inplace=True)

        # combine features
        df_feature = pd.concat(l_df_feature)
        df_feature.drop_duplicates(inplace=True)

        Write_10X(df_mtx, df_feature, path_folder_mtx_10x_output)

        # split a matrix file into multiple files
        MTX_10X_Split(
            path_folder_mtx_10x_output,
            int_max_num_entries_for_chunk=int_max_num_entries_for_chunk,
        )


def __Combine_Dictionaries__(path_folder_mtx_10x_input, name_dict):
    """# 2022-03-06 00:06:23
    combined dictionaries processed from individual files
    """
    import collections

    if filesystem_operations(
        "exists", f"{path_folder_mtx_10x_input}{name_dict}.tsv.gz"
    ):
        """if an output file already exists, read the file and return the combined dictionary"""
        dict_combined = (
            pd.read_csv(
                f"{path_folder_mtx_10x_input}{name_dict}.tsv.gz",
                sep="\t",
                header=None,
                index_col=0,
            )
            .iloc[:, 0]
            .to_dict()
        )
    else:
        """combine summarized results"""
        l_path_file = filesystem_operations(
            "glob", f"{path_folder_mtx_10x_input}{name_dict}.*"
        )
        try:
            counter = collections.Counter(
                pd.read_csv(l_path_file[0], sep="\t", header=None, index_col=0)
                .iloc[:, 0]
                .to_dict()
            )  # initialize counter object with the dictionary from the first file
        except pd.errors.EmptyDataError:
            counter = (
                collections.Counter()
            )  # when an error (possibly because the file is empty) occur, use an empty counter
        for path_file in l_path_file[1:]:
            # when an error (possibly because the file is empty) occur, skip updating the counter
            try:
                counter = counter + collections.Counter(
                    pd.read_csv(path_file, sep="\t", header=None, index_col=0)
                    .iloc[:, 0]
                    .to_dict()
                )  # update counter object using the dictionary from each file
            except pd.errors.EmptyDataError:
                pass
        dict_combined = dict(counter)  # retrieve a combined dictionary
        """remove temporary files """
        for path_file in l_path_file:
            filesystem_operations("rm", path_file)
        """ save dictionary as a file """
        pd.Series(dict_combined).to_csv(
            f"{path_folder_mtx_10x_input}{name_dict}.tsv.gz", sep="\t", header=None
        )
    return dict_combined  # returns a combined dictionary


def __MTX_10X_Summarize_Counts__summarize_counts_for_each_mtx_10x__(
    path_file_input: str, path_folder_mtx_10x_input: str, str_data_type: str
):
    """
    internal function for MTX_10X_Summarize_Count
    # 2023-09-10 16:37:13

    str_data_type : str # matrix market datatype in string format
    """
    """ survey the metrics """
    """ for each split mtx file, count number of umi and n_feature for each cells or the number of cells for each feature """
    """ initialize the dictionaries that will be handled by the current function """
    dict_id_column_to_count = dict()
    dict_id_column_to_n_features = dict()
    dict_id_row_to_count = dict()
    dict_id_row_to_n_cells = dict()
    dict_id_row_to_log_transformed_count = dict()

    flag_matrix_contain_float_values = (
        str_data_type == "real"
    )  # retrieve a flag indicating the matrix is containing float values

    global dict_name_set_feature_to_set_id_row  # use global read-only object
    dict_name_set_feature_to_dict_id_column_to_count = dict(
        (name_set_feature, dict())
        for name_set_feature in dict_name_set_feature_to_set_id_row
    )  # initialize 'dict_name_set_feature_to_dict_id_column_to_count'
    for path_file_input_mtx in pd.read_csv(
        path_file_input, sep="\t", header=None
    ).values.ravel():
        with gzip.open(path_file_input_mtx, "rb") as file:
            """read the first line"""
            line = file.readline().decode()
            """ if the first line of the file contains a comment line, read all comment lines and a description line following the comments. """
            if len(line) > 0 and line[0] == "%":
                # read comment and the description line
                while True:
                    if line[0] != "%":
                        break
                    line = file.readline().decode()  # read the next line
                # process the description line
                int_num_rows, int_num_columns, int_num_entries = tuple(
                    int(e) for e in line.strip().split()
                )  # retrieve the number of rows, number of columns and number of entries
                line = file.readline().decode()  # read the next line
            """ process entries"""
            while True:
                if len(line) == 0:
                    break
                """ parse a record, and update metrics """
                l_str = (
                    line.strip().split()
                )  # parse a record of a matrix-market format file
                if flag_matrix_contain_float_values:  # parse float data type
                    id_row, id_column, value = (
                        int(l_str[0]),
                        int(l_str[1]),
                        float(l_str[2]),
                    )
                else:  # parse integer data type
                    id_row, id_column, value = tuple(int(float(e)) for e in l_str)

                """ 1-based > 0-based coordinates """
                id_row -= 1
                id_column -= 1
                """ update umi count for each cell """
                if id_column not in dict_id_column_to_count:
                    dict_id_column_to_count[id_column] = 0
                dict_id_column_to_count[id_column] += value
                """ update umi count of specific set of features for each cell """
                for (
                    name_set_feature
                ) in dict_name_set_feature_to_dict_id_column_to_count:
                    if id_row in dict_name_set_feature_to_set_id_row[name_set_feature]:
                        if (
                            id_column
                            not in dict_name_set_feature_to_dict_id_column_to_count[
                                name_set_feature
                            ]
                        ):
                            dict_name_set_feature_to_dict_id_column_to_count[
                                name_set_feature
                            ][id_column] = 0
                        dict_name_set_feature_to_dict_id_column_to_count[
                            name_set_feature
                        ][id_column] += value
                """ update n_features for each cell """
                if id_column not in dict_id_column_to_n_features:
                    dict_id_column_to_n_features[id_column] = 0
                dict_id_column_to_n_features[id_column] += 1
                """ update umi count for each feature """
                if id_row not in dict_id_row_to_count:
                    dict_id_row_to_count[id_row] = 0
                dict_id_row_to_count[id_row] += value
                """ update n_cells for each feature """
                if id_row not in dict_id_row_to_n_cells:
                    dict_id_row_to_n_cells[id_row] = 0
                dict_id_row_to_n_cells[id_row] += 1
                """ update log transformed counts, calculated by 'X_new = log_10(X_old + 1)', for each feature """
                if id_row not in dict_id_row_to_log_transformed_count:
                    dict_id_row_to_log_transformed_count[id_row] = 0
                dict_id_row_to_log_transformed_count[id_row] += math.log10(value + 1)

                """ read the next line """
                line = (
                    file.readline().decode()
                )  # binary > uncompressed string # read the next line

    # save collected count as tsv files
    str_uuid_process = bk.UUID()  # retrieve uuid of the current process
    pd.Series(dict_id_column_to_count).to_csv(
        f"{path_folder_mtx_10x_input}dict_id_column_to_count.{str_uuid_process}.tsv.gz",
        sep="\t",
        header=None,
    )
    pd.Series(dict_id_column_to_n_features).to_csv(
        f"{path_folder_mtx_10x_input}dict_id_column_to_n_features.{str_uuid_process}.tsv.gz",
        sep="\t",
        header=None,
    )
    pd.Series(dict_id_row_to_count).to_csv(
        f"{path_folder_mtx_10x_input}dict_id_row_to_count.{str_uuid_process}.tsv.gz",
        sep="\t",
        header=None,
    )
    pd.Series(dict_id_row_to_n_cells).to_csv(
        f"{path_folder_mtx_10x_input}dict_id_row_to_n_cells.{str_uuid_process}.tsv.gz",
        sep="\t",
        header=None,
    )
    pd.Series(dict_id_row_to_log_transformed_count).to_csv(
        f"{path_folder_mtx_10x_input}dict_id_row_to_log_transformed_count.{str_uuid_process}.tsv.gz",
        sep="\t",
        header=None,
    )

    # save collected counts as tsv files for 'dict_name_set_feature_to_dict_id_column_to_count'
    for name_set_feature in dict_name_set_feature_to_dict_id_column_to_count:
        pd.Series(
            dict_name_set_feature_to_dict_id_column_to_count[name_set_feature]
        ).to_csv(
            f"{path_folder_mtx_10x_input}{name_set_feature}.dict_id_column_to_count.{str_uuid_process}.tsv.gz",
            sep="\t",
            header=None,
        )


def MTX_10X_Summarize_Counts(
    path_folder_mtx_10x_input,
    verbose=False,
    int_num_threads=15,
    flag_split_mtx=True,
    int_max_num_entries_for_chunk=10000000,
    dict_name_set_feature_to_l_id_feature=dict(),
    flag_split_mtx_again=False,
):
    """# 2022-04-28 06:53:45
    Summarize
    (1) UMI and Feature counts for each cell,
    (2) UMI and Cell counts for each feature, and
    (3) log10-transformed values of UMI counts (X_new = log_10(X_old + 1)) for each feature
    (4) UMI counts for the optionally given lists of features for each cell
    and save these metrics as TSV files

    Inputs:
    'dict_name_set_feature_to_l_id_feature' : (Default: None)
                                            a dictionary with 'name_set_features' as key and a list of id_feature as value for each set of id_features.
                                            If None is given, only the basic metrics will be summarized.
                                            'name_set_features' should be compatible as a Linux file system (should not contain '/' and other special characters, such as newlines).
                                            (for Scarab short_read outputs)
                                            If 'atac' is given, 'promoter_and_gene_body', 'promoter' features will be summarized.
                                            If 'multiome' is given, total 'atac' counts will be summarized separately in addition to 'atac' mode

    Returns:
    a dictionary containing the following and other additional dictionaries: dict_id_column_to_count, dict_id_column_to_n_features, dict_id_row_to_count, dict_id_row_to_n_cells, dict_id_row_to_log_transformed_count
    """

    """ the name of the dictionaries handled by this function (basic) """
    l_name_dict = [
        "dict_id_column_to_count",
        "dict_id_column_to_n_features",
        "dict_id_row_to_count",
        "dict_id_row_to_n_cells",
        "dict_id_row_to_log_transformed_count",
    ]

    """ handle inputs """
    if path_folder_mtx_10x_input[-1] != "/":
        path_folder_mtx_10x_input += "/"

    # define flag and check whether the flag exists
    path_file_flag = f"{path_folder_mtx_10x_input}counts_summarized.flag"
    if not filesystem_operations("exists", path_file_flag):
        # define input file directories
        path_file_input_bc = f"{path_folder_mtx_10x_input}barcodes.tsv.gz"
        path_file_input_feature = f"{path_folder_mtx_10x_input}features.tsv.gz"
        path_file_input_mtx = f"{path_folder_mtx_10x_input}matrix.mtx.gz"

        # check whether all required files are present
        if sum(
            list(
                not filesystem_operations("exists", path_folder)
                for path_folder in [
                    path_file_input_bc,
                    path_file_input_feature,
                    path_file_input_mtx,
                ]
            )
        ):
            if verbose:
                logger.info(f"required file(s) is not present at {path_folder_mtx_10x}")

        """ split input mtx file into multiple files """
        l_path_file_mtx_10x = MTX_10X_Split(
            path_folder_mtx_10x_input,
            int_max_num_entries_for_chunk=int_max_num_entries_for_chunk,
            flag_split_mtx=flag_split_mtx,
            flag_split_mtx_again=flag_split_mtx_again,
        )
        str_datatype = _MTX_Detect_data_type(
            path_file_input_mtx
        )  # retrieve the data type of the matrix in string format

        """ prepare 'dict_name_set_feature_to_set_id_row' for summarizing total counts for given sets of features """
        global dict_name_set_feature_to_set_id_row
        dict_name_set_feature_to_set_id_row = (
            dict()
        )  # initialize 'dict_name_set_feature_to_set_id_row'
        if dict_name_set_feature_to_l_id_feature is not None:
            arr_id_feature = pd.read_csv(
                path_file_input_feature, sep="\t", usecols=[0], header=None
            ).values.ravel()  # retrieve array of id_features
            dict_id_feature_to_id_row = dict(
                (e, i) for i, e in enumerate(arr_id_feature)
            )  # retrieve id_feature -> id_row mapping

            """ handle presets for 'dict_name_set_feature_to_l_id_feature' """
            if isinstance(dict_name_set_feature_to_l_id_feature, str):
                str_preset = dict_name_set_feature_to_l_id_feature  # retrieve preset
                dict_name_set_feature_to_l_id_feature = (
                    dict()
                )  # initialize the dictionary
                if str_preset in ["multiome", "atac"]:
                    if str_preset == "multiome":
                        arr_id_feature_atac = (
                            bk.Search_list_of_strings_with_multiple_query(
                                arr_id_feature, "|mode=atac"
                            )
                        )
                        dict_name_set_feature_to_l_id_feature["atac_all"] = (
                            arr_id_feature_atac
                        )
                    elif str_preset == "atac":
                        arr_id_feature_atac = arr_id_feature
                    # add sets of promoter and gene_body features
                    arr_id_feature_atac_promoter_and_gene_body = (
                        bk.Search_list_of_strings_with_multiple_query(
                            arr_id_feature_atac,
                            "-genomic_region|",
                            "-repeatmasker_ucsc|",
                            "-regulatory_element|",
                        )
                    )
                    arr_id_feature_atac_promoter = (
                        bk.Search_list_of_strings_with_multiple_query(
                            arr_id_feature_atac_promoter_and_gene_body, "promoter|"
                        )
                    )
                    dict_name_set_feature_to_l_id_feature[
                        "atac_promoter_and_gene_body"
                    ] = arr_id_feature_atac_promoter_and_gene_body
                    dict_name_set_feature_to_l_id_feature["atac_promoter"] = (
                        arr_id_feature_atac_promoter
                    )

            # make sure that 'name_set_feature' does not contains characters incompatible with linux file path
            for name_set_feature in dict_name_set_feature_to_l_id_feature:
                assert not ("/" in name_set_feature or "\n" in name_set_feature)

            dict_name_set_feature_to_set_id_row = dict(
                (
                    name_set_feature,
                    set(
                        dict_id_feature_to_id_row[id_feature]
                        for id_feature in dict_name_set_feature_to_l_id_feature[
                            name_set_feature
                        ]
                    ),
                )
                for name_set_feature in dict_name_set_feature_to_l_id_feature
            )
            # bk.PICKLE_Write( f"{path_folder_mtx_10x_input}dict_name_set_feature_to_set_id_row.binary.pickle", dict_name_set_feature_to_set_id_row ) # write the dictionary as a pickle

        """ summarize each split mtx file """
        bk.Multiprocessing(
            l_path_file_mtx_10x,
            __MTX_10X_Summarize_Counts__summarize_counts_for_each_mtx_10x__,
            n_threads=int_num_threads,
            global_arguments=[path_folder_mtx_10x_input, str_datatype],
        )

        """ combine summarized results """
        # update the list of the names of dictionaries
        l_name_dict += list(
            f"{name_set_feature}.dict_id_column_to_count"
            for name_set_feature in bk.GLOB_Retrive_Strings_in_Wildcards(
                f"{path_folder_mtx_10x_input}*.dict_id_column_to_count.*.tsv.gz"
            ).wildcard_0.unique()
        )

        dict_dict = dict()
        for name_dict in l_name_dict:
            dict_dict[name_dict] = __Combine_Dictionaries__(
                path_folder_mtx_10x_input, name_dict
            )
        # write the flag
        with open(path_file_flag, "w") as newfile:
            newfile.write("completed at " + bk.TIME_GET_timestamp(True))
    else:
        """read summarized results"""
        # update the list of the names of dictionaries
        l_name_dict += list(
            f"{name_set_feature}.dict_id_column_to_count"
            for name_set_feature in bk.GLOB_Retrive_Strings_in_Wildcards(
                f"{path_folder_mtx_10x_input}*.dict_id_column_to_count.tsv.gz"
            ).wildcard_0.unique()
        )

        dict_dict = dict()
        for name_dict in l_name_dict:
            try:
                dict_dict[name_dict] = (
                    pd.read_csv(
                        f"{path_folder_mtx_10x_input}{name_dict}.tsv.gz",
                        sep="\t",
                        header=None,
                        index_col=0,
                    )
                    .iloc[:, 0]
                    .to_dict()
                )
            except (
                pd.errors.EmptyDataError
            ):  # handle when the current dictionary is empty
                dict_dict[name_dict] = dict()

    # return summarized metrics
    return dict_dict


def MTX_10X_Retrieve_number_of_rows_columns_and_records(path_folder_mtx_10x_input):
    """# 2022-03-05 19:58:32
    Retrieve the number of rows, columns, and entries from the matrix with the matrix market format

    'path_folder_mtx_10x_input' : a folder mtx file resides or path to mtx file

    Returns:
    int_num_rows, int_num_columns, int_num_entries
    """
    """ handle inputs """
    if (
        path_folder_mtx_10x_input[-3:].lower() == ".gz"
    ):  # when a path to mtx file was given
        path_file_input_mtx = path_folder_mtx_10x_input
    else:  # when a folder where mtx file resides was given
        if path_folder_mtx_10x_input[-1] != "/":
            path_folder_mtx_10x_input += "/"

        # define input file directories
        path_file_input_mtx = f"{path_folder_mtx_10x_input}matrix.mtx.gz"

        # check whether all required files are present
        if sum(
            list(
                not filesystem_operations("exists", path_folder)
                for path_folder in [path_file_input_mtx]
            )
        ):
            return None

    # read the input matrix
    with gzip.open(path_file_input_mtx, "rb") as file:
        """read the first line"""
        line = file.readline().decode().strip()
        """ if the first line of the file contains a comment line, read all comment lines and a description line following the comments. """
        if len(line) > 0 and line[0] == "%":
            # read comment and the description line
            while True:
                if line[0] != "%":
                    break
                line = file.readline().decode().strip()  # read the next line
            # process the description line
            int_num_rows, int_num_columns, int_num_entries = tuple(
                int(e) for e in line.strip().split()
            )  # retrieve the number of rows, number of columns and number of entries
        else:
            """the first line does not contain a comment, assumes it contains a description line"""
            int_num_rows, int_num_columns, int_num_entries = tuple(
                int(e) for e in line.strip().split()
            )  # retrieve the number of rows, number of columns and number of entries
    return int_num_rows, int_num_columns, int_num_entries


(
    dict_id_column_to_count,
    dict_id_row_to_avg_count,
    dict_id_row_to_avg_log_transformed_count,
    dict_id_row_to_avg_normalized_count,
    dict_id_row_to_avg_log_transformed_normalized_count,
) = (
    dict(),
    dict(),
    dict(),
    dict(),
    dict(),
)  # global variables # total UMI counts for each cell, average feature counts for each feature


def __MTX_10X_Calculate_Average_Log10_Transformed_Normalized_Expr__first_pass__(
    path_file_input, path_folder_mtx_10x_input, int_target_sum
):
    """# 2022-03-06 01:21:07
    internal function for MTX_10X_Calculate_Average_Log10_Transformed_Normalized_Expr
    """
    global dict_id_column_to_count, dict_id_row_to_avg_count, dict_id_row_to_avg_log_transformed_count  # use data in read-only global variables
    """ initialize the dictionaries that will be handled by the current function """
    dict_id_row_to_deviation_from_mean_count = dict()
    dict_id_row_to_deviation_from_mean_log_transformed_count = dict()
    dict_id_row_to_normalized_count = dict()
    dict_id_row_to_log_transformed_normalized_count = dict()

    for path_file_input_mtx in pd.read_csv(
        path_file_input, sep="\t", header=None
    ).values.ravel():
        with gzip.open(path_file_input_mtx, "rb") as file:
            """read the first line"""
            line = file.readline().decode()
            """ if the first line of the file contains a comment line, read all comment lines and a description line following the comments. """
            if len(line) > 0 and line[0] == "%":
                # read comment and the description line
                while True:
                    if line[0] != "%":
                        break
                    line = file.readline().decode()  # read the next line
                # process the description line
                int_num_rows, int_num_columns, int_num_entries = tuple(
                    int(e) for e in line.strip().split()
                )  # retrieve the number of rows, number of columns and number of entries
                line = file.readline().decode()  # read the next line
            """ process entries"""
            while True:
                if len(line) == 0:
                    break
                """ parse a record, and update metrics """
                id_row, id_column, int_value = tuple(
                    int(e) for e in line.strip().split()
                )  # parse a record of a matrix-market format file
                """ 1-based > 0-based coordinates """
                id_row -= 1
                id_column -= 1

                """ update deviation from mean umi count for count of each feature """
                if id_row not in dict_id_row_to_deviation_from_mean_count:
                    dict_id_row_to_deviation_from_mean_count[id_row] = 0
                dict_id_row_to_deviation_from_mean_count[id_row] += (
                    int_value - dict_id_row_to_avg_count[id_row]
                ) ** 2
                """ update deviation from mean log transformed umi count for log_transformed count of each feature """
                if (
                    id_row
                    not in dict_id_row_to_deviation_from_mean_log_transformed_count
                ):
                    dict_id_row_to_deviation_from_mean_log_transformed_count[id_row] = 0
                dict_id_row_to_deviation_from_mean_log_transformed_count[id_row] += (
                    math.log10(int_value + 1)
                    - dict_id_row_to_avg_log_transformed_count[id_row]
                ) ** 2
                """ calculate normalized target sum """
                int_value_normalized = (
                    int_value / dict_id_column_to_count[id_column] * int_target_sum
                )
                """ update normalized counts, calculated by 'X_new = X_old / total_umi * int_target_sum', for each feature """
                if id_row not in dict_id_row_to_normalized_count:
                    dict_id_row_to_normalized_count[id_row] = 0
                dict_id_row_to_normalized_count[id_row] += int_value_normalized
                """ update log transformed normalized counts, calculated by 'X_new = log_10(X_old / total_umi * int_target_sum + 1)', for each feature """
                if id_row not in dict_id_row_to_log_transformed_normalized_count:
                    dict_id_row_to_log_transformed_normalized_count[id_row] = 0
                dict_id_row_to_log_transformed_normalized_count[id_row] += math.log10(
                    int_value_normalized + 1
                )

                line = (
                    file.readline().decode()
                )  # binary > uncompressed string # read the next line

    # save collected count as tsv files
    str_uuid_process = bk.UUID()  # retrieve uuid of the current process
    pd.Series(dict_id_row_to_deviation_from_mean_count).to_csv(
        f"{path_folder_mtx_10x_input}dict_id_row_to_deviation_from_mean_count.{str_uuid_process}.tsv.gz",
        sep="\t",
        header=None,
    )
    pd.Series(dict_id_row_to_deviation_from_mean_log_transformed_count).to_csv(
        f"{path_folder_mtx_10x_input}dict_id_row_to_deviation_from_mean_log_transformed_count.{str_uuid_process}.tsv.gz",
        sep="\t",
        header=None,
    )
    pd.Series(dict_id_row_to_normalized_count).to_csv(
        f"{path_folder_mtx_10x_input}dict_id_row_to_normalized_count.{str_uuid_process}.tsv.gz",
        sep="\t",
        header=None,
    )
    pd.Series(dict_id_row_to_log_transformed_normalized_count).to_csv(
        f"{path_folder_mtx_10x_input}dict_id_row_to_log_transformed_normalized_count.{str_uuid_process}.tsv.gz",
        sep="\t",
        header=None,
    )


def __MTX_10X_Calculate_Average_Log10_Transformed_Normalized_Expr__second_pass__(
    path_file_input, path_folder_mtx_10x_input, int_target_sum
):
    """# 2022-03-06 01:21:14
    internal function for MTX_10X_Calculate_Average_Log10_Transformed_Normalized_Expr
    """
    global dict_id_column_to_count, dict_id_row_to_avg_normalized_count, dict_id_row_to_avg_log_transformed_normalized_count  # use data in read-only global variables
    """ initialize the dictionaries that will be handled by the current function """
    dict_id_row_to_deviation_from_mean_normalized_count = dict()
    dict_id_row_to_deviation_from_mean_log_transformed_normalized_count = dict()

    for path_file_input_mtx in pd.read_csv(
        path_file_input, sep="\t", header=None
    ).values.ravel():
        with gzip.open(path_file_input_mtx, "rb") as file:
            """read the first line"""
            line = file.readline().decode()
            """ if the first line of the file contains a comment line, read all comment lines and a description line following the comments. """
            if len(line) > 0 and line[0] == "%":
                # read comment and the description line
                while True:
                    if line[0] != "%":
                        break
                    line = file.readline().decode()  # read the next line
                # process the description line
                int_num_rows, int_num_columns, int_num_entries = tuple(
                    int(e) for e in line.strip().split()
                )  # retrieve the number of rows, number of columns and number of entries
                line = file.readline().decode()  # read the next line
            """ process entries"""
            while True:
                if len(line) == 0:
                    break
                """ parse a record, and update metrics """
                id_row, id_column, int_value = tuple(
                    int(e) for e in line.strip().split()
                )  # parse a record of a matrix-market format file
                """ 1-based > 0-based coordinates """
                id_row -= 1
                id_column -= 1

                """ calculate normalized target sum """
                int_value_normalized = (
                    int_value / dict_id_column_to_count[id_column] * int_target_sum
                )
                """ update deviation from mean normalized umi counts, calculated by 'X_new = X_old / total_umi * int_target_sum', for each feature """
                if id_row not in dict_id_row_to_deviation_from_mean_normalized_count:
                    dict_id_row_to_deviation_from_mean_normalized_count[id_row] = 0
                dict_id_row_to_deviation_from_mean_normalized_count[id_row] += (
                    int_value_normalized - dict_id_row_to_avg_normalized_count[id_row]
                ) ** 2
                """ update deviation from mean log transformed normalized umi counts, calculated by 'X_new = log_10(X_old / total_umi * int_target_sum + 1)', for each feature """
                if (
                    id_row
                    not in dict_id_row_to_deviation_from_mean_log_transformed_normalized_count
                ):
                    dict_id_row_to_deviation_from_mean_log_transformed_normalized_count[
                        id_row
                    ] = 0
                dict_id_row_to_deviation_from_mean_log_transformed_normalized_count[
                    id_row
                ] += (
                    math.log10(int_value_normalized + 1)
                    - dict_id_row_to_avg_log_transformed_normalized_count[id_row]
                ) ** 2

                line = (
                    file.readline().decode()
                )  # binary > uncompressed string # read the next line

    # save collected count as tsv files
    str_uuid_process = bk.UUID()  # retrieve uuid of the current process
    pd.Series(dict_id_row_to_deviation_from_mean_normalized_count).to_csv(
        f"{path_folder_mtx_10x_input}dict_id_row_to_deviation_from_mean_normalized_count.{str_uuid_process}.tsv.gz",
        sep="\t",
        header=None,
    )
    pd.Series(
        dict_id_row_to_deviation_from_mean_log_transformed_normalized_count
    ).to_csv(
        f"{path_folder_mtx_10x_input}dict_id_row_to_deviation_from_mean_log_transformed_normalized_count.{str_uuid_process}.tsv.gz",
        sep="\t",
        header=None,
    )


def MTX_10X_Calculate_Average_Log10_Transformed_Normalized_Expr(
    path_folder_mtx_10x_input,
    int_target_sum=10000,
    verbose=False,
    int_num_threads=15,
    flag_split_mtx=True,
    int_max_num_entries_for_chunk=10000000,
):
    """# 2022-02-23 22:54:35
    Calculate average log transformed normalized expr
    (1) UMI and Feature counts for cells, and
    (2) Cell counts for features,
    and save these metrics as TSV files

    Arguments:
    'int_target_sum' : the target count for the total UMI count for each cell. The counts will normalized to meet the target sum.

    Returns:
    dict_id_column_to_count, dict_id_column_to_n_features, dict_id_row_to_count, dict_id_row_to_n_cells, dict_id_row_to_log_transformed_count
    """

    """ handle inputs """
    if path_folder_mtx_10x_input[-1] != "/":
        path_folder_mtx_10x_input += "/"

    # define flag and check whether the flag exists
    path_file_flag = f"{path_folder_mtx_10x_input}avg_expr_normalized_summarized.int_target_sum__{int_target_sum}.flag"
    if not filesystem_operations("exists", path_file_flag):
        # define input file directories
        path_file_input_bc = f"{path_folder_mtx_10x_input}barcodes.tsv.gz"
        path_file_input_feature = f"{path_folder_mtx_10x_input}features.tsv.gz"
        path_file_input_mtx = f"{path_folder_mtx_10x_input}matrix.mtx.gz"

        # check whether all required files are present
        if sum(
            list(
                not filesystem_operations("exists", path_folder)
                for path_folder in [
                    path_file_input_bc,
                    path_file_input_feature,
                    path_file_input_mtx,
                ]
            )
        ):
            if verbose:
                logger.info(f"required file(s) is not present at {path_folder_mtx_10x}")

        """ split input mtx file into multiple files """
        l_path_file_mtx_10x = MTX_10X_Split(
            path_folder_mtx_10x_input,
            int_max_num_entries_for_chunk=int_max_num_entries_for_chunk,
            flag_split_mtx=flag_split_mtx,
        )

        """ retrieve number of cells, features, and entries from the matrix file """
        (
            int_num_cells,
            int_num_features,
            int_num_entries,
        ) = MTX_10X_Retrieve_number_of_rows_columns_and_records(
            path_folder_mtx_10x_input
        )

        """ summarizes counts """
        global dict_id_column_to_count, dict_id_row_to_avg_count, dict_id_row_to_avg_log_transformed_count, dict_id_row_to_avg_normalized_count, dict_id_row_to_avg_log_transformed_normalized_count  # use global variable
        dict_data = MTX_10X_Summarize_Counts(
            path_folder_mtx_10x_input,
            verbose=verbose,
            int_num_threads=int_num_threads,
            flag_split_mtx=flag_split_mtx,
            int_max_num_entries_for_chunk=int_max_num_entries_for_chunk,
        )
        (
            dict_id_column_to_count,
            dict_id_column_to_n_features,
            dict_id_row_to_count,
            dict_id_row_to_n_cells,
            dict_id_row_to_log_transformed_count,
        ) = (
            dict_data["dict_id_column_to_count"],
            dict_data["dict_id_column_to_n_features"],
            dict_data["dict_id_row_to_count"],
            dict_data["dict_id_row_to_n_cells"],
            dict_data["dict_id_row_to_log_transformed_count"],
        )  # parse 'dict_data'

        """ first pass """
        # calculate mean counts
        dict_id_row_to_avg_count = (
            pd.Series(dict_id_row_to_count) / int_num_cells
        ).to_dict()  # calculate average expression of each feature
        dict_id_row_to_avg_log_transformed_count = (
            pd.Series(dict_id_row_to_log_transformed_count) / int_num_cells
        ).to_dict()  # calculate average log-transformed expression of each feature

        """ calculated average log2 transformed normalized expr for each split mtx file """
        bk.Multiprocessing(
            l_path_file_mtx_10x,
            __MTX_10X_Calculate_Average_Log10_Transformed_Normalized_Expr__first_pass__,
            n_threads=int_num_threads,
            global_arguments=[path_folder_mtx_10x_input, int_target_sum],
        )

        l_name_dict_first_pass = [
            "dict_id_row_to_deviation_from_mean_count",
            "dict_id_row_to_deviation_from_mean_log_transformed_count",
            "dict_id_row_to_normalized_count",
            "dict_id_row_to_log_transformed_normalized_count",
        ]

        """ combine summarized results """
        dict_dict = dict()
        for name_dict in l_name_dict_first_pass:
            dict_dict[name_dict] = __Combine_Dictionaries__(
                path_folder_mtx_10x_input, name_dict
            )

        """ second pass """
        # calculate mean counts
        dict_id_row_to_avg_normalized_count = (
            pd.Series(dict_dict["dict_id_row_to_normalized_count"]) / int_num_cells
        ).to_dict()  # calculate average expression of each feature
        dict_id_row_to_avg_log_transformed_normalized_count = (
            pd.Series(dict_dict["dict_id_row_to_log_transformed_normalized_count"])
            / int_num_cells
        ).to_dict()  # calculate average log-transformed expression of each feature

        """ calculated average log2 transformed normalized expr for each split mtx file """
        bk.Multiprocessing(
            l_path_file_mtx_10x,
            __MTX_10X_Calculate_Average_Log10_Transformed_Normalized_Expr__second_pass__,
            n_threads=int_num_threads,
            global_arguments=[path_folder_mtx_10x_input, int_target_sum],
        )

        l_name_dict_second_pass = [
            "dict_id_row_to_deviation_from_mean_normalized_count",
            "dict_id_row_to_deviation_from_mean_log_transformed_normalized_count",
        ]

        """ combine summarized results """
        for name_dict in l_name_dict_second_pass:
            dict_dict[name_dict] = __Combine_Dictionaries__(
                path_folder_mtx_10x_input, name_dict
            )

        """ compose a dataframe containing the summary about the features """
        df_summary = pd.DataFrame(
            {
                "n_cells": pd.Series(dict_id_row_to_n_cells),
                "variance_of_count": pd.Series(
                    dict_dict["dict_id_row_to_deviation_from_mean_count"]
                )
                / (int_num_cells - 1),
                "variance_of_log_transformed_count": pd.Series(
                    dict_dict[
                        "dict_id_row_to_deviation_from_mean_log_transformed_count"
                    ]
                )
                / (int_num_cells - 1),
                "variance_of_normalized_count": pd.Series(
                    dict_dict["dict_id_row_to_deviation_from_mean_normalized_count"]
                )
                / (int_num_cells - 1),
                "variance_of_log_transformed_normalized_count": pd.Series(
                    dict_dict[
                        "dict_id_row_to_deviation_from_mean_log_transformed_normalized_count"
                    ]
                )
                / (int_num_cells - 1),
                "mean_count": pd.Series(dict_id_row_to_avg_count),
                "mean_log_transformed_count": pd.Series(
                    dict_id_row_to_avg_log_transformed_count
                ),
                "mean_normalized_count": pd.Series(dict_id_row_to_avg_normalized_count),
                "mean_log_transformed_normalized_count": pd.Series(
                    dict_id_row_to_avg_log_transformed_normalized_count
                ),
            }
        )
        # read a dataframe containing features
        df_feature = pd.read_csv(path_file_input_feature, sep="\t", header=None)
        df_feature.columns = ["id_feature", "feature", "feature_type"]

        df_summary = df_summary.join(
            df_feature, how="left"
        )  # add df_feature to the df_summary
        df_summary.index.name = "id_row"
        df_summary.reset_index(drop=False, inplace=True)  # retrieve id_row as a column
        df_summary.to_csv(
            f"{path_folder_mtx_10x_input}statistical_summary_of_features.int_target_sum__{int_target_sum}.tsv.gz",
            sep="\t",
            index=False,
        )  # save statistical summary as a text file

        # write the flag
        with open(path_file_flag, "w") as newfile:
            newfile.write("completed at " + bk.TIME_GET_timestamp(True))
    else:
        """if 'MTX_10X_Calculate_Average_Log10_Transformed_Normalized_Expr' function has been already run on the current folder, read the previously saved result, and return the summary dataframe"""
        df_summary = pd.read_csv(
            f"{path_folder_mtx_10x_input}statistical_summary_of_features.int_target_sum__{int_target_sum}.tsv.gz",
            sep="\t",
        )  # save statistical summary as a text file
    return df_summary


dict_id_column_previous_to_id_column_current, dict_id_row_previous_to_id_row_current = (
    dict(),
    dict(),
)


def __MTX_10X_Filter__filter_mtx_10x__(
    path_file_input: str,
    path_folder_mtx_10x_output: str,
    str_data_type: str,
):
    """# 2023-09-10 16:37:33
    __MTX_10X_Filter__filter_mtx_10x__

    str_data_type : str # matrix market datatype in string format

    Returns:
    int_n_entries = total number of entries written by the current process after filtering
    """
    int_n_entries = (
        0  # total number of entries written by the current process after filtering
    )
    flag_matrix_contain_float_values = (
        str_data_type == "real"
    )  # retrieve a flag indicating the matrix is containing float values
    #     dict_id_column_previous_to_id_column_current = bk.PICKLE_Read( f'{path_folder_mtx_10x_output}dict_id_column_previous_to_id_column_current.pickle' ) # retrieve id_feature to index_feature mapping
    #     dict_id_row_previous_to_id_row_current = bk.PICKLE_Read( f'{path_folder_mtx_10x_output}dict_id_row_previous_to_id_row_current.pickle' ) # retrieve id_feature to index_feature mapping
    """ write a filtered matrix.mtx.gz for each split mtx file """
    for path_file_mtx_10x, index_mtx_10x in pd.read_csv(
        path_file_input, sep="\t"
    ).values:
        # directly write matrix.mtx.gz file without using an external dependency
        with gzip.open(
            f"{path_folder_mtx_10x_output}matrix.mtx.gz.{index_mtx_10x}.gz", "wb"
        ) as newfile:
            with gzip.open(path_file_mtx_10x, "rb") as file:
                """read the first line"""
                line = file.readline().decode()
                """ if the first line of the file contains a comment line, read all comment lines and a description line following the comments. """
                if len(line) > 0 and line[0] == "%":
                    # read comment and the description line
                    while True:
                        if line[0] != "%":
                            break
                        line = file.readline().decode()  # read the next line
                    # process the description line
                    int_num_rows, int_num_columns, int_num_entries = tuple(
                        int(e) for e in line.strip().split()
                    )  # retrieve the number of rows, number of columns and number of entries
                    line = file.readline().decode()  # read the next line
                """ process entries"""
                while True:
                    if len(line) == 0:
                        break
                    """ parse a record """
                    l_str = (
                        line.strip().split()
                    )  # parse a record of a matrix-market format file
                    if flag_matrix_contain_float_values:  # parse float data type
                        id_row, id_column, value = (
                            int(l_str[0]),
                            int(l_str[1]),
                            float(l_str[2]),
                        )
                    else:  # parse integer data type
                        id_row, id_column, value = tuple(int(float(e)) for e in l_str)

                    """ 1-based > 0-based coordinates """
                    id_row -= 1
                    id_column -= 1
                    """ write a record to the new matrix file only when both id_row and id_column belongs to filtered id_rows and id_columns """
                    if (
                        id_row in dict_id_row_previous_to_id_row_current
                        and id_column in dict_id_column_previous_to_id_column_current
                    ):
                        newfile.write(
                            (
                                " ".join(
                                    tuple(
                                        map(
                                            str,
                                            [
                                                dict_id_row_previous_to_id_row_current[
                                                    id_row
                                                ]
                                                + 1,
                                                dict_id_column_previous_to_id_column_current[
                                                    id_column
                                                ]
                                                + 1,
                                                value,
                                            ],
                                        )
                                    )
                                )
                                + "\n"
                            ).encode()
                        )  # map id_row and id_column of the previous matrix to those of the filtered matrix (new matrix) # 0-based > 1-based coordinates
                        int_n_entries += 1  # update the total number of entries written by the current process
                    line = file.readline().decode()  # read the next line
    return int_n_entries  # returns the total number of entries written by the current process


def MTX_10X_Filter(
    path_folder_mtx_10x_input,
    path_folder_mtx_10x_output,
    min_counts=None,
    min_features=None,
    min_cells=None,
    l_features=None,
    l_cells=None,
    verbose=False,
    function_for_adjusting_thresholds=None,
    int_num_threads=15,
    flag_split_mtx=True,
    int_max_num_entries_for_chunk=10000000,
):
    """# 2022-08-20 10:23:28
    # hyunsu-an
    read 10x count matrix and filter matrix based on several thresholds
    'path_folder_mtx_10x_input' : a folder containing files for the input 10x count matrix
    'path_folder_mtx_10x_output' : a folder containing files for the input 10x count matrix

    Only the threshold arguments for either cells ( 'min_counts', 'min_features' ) or features ( 'min_cells' ) can be given at a time.

    'min_counts' : the minimum number of total counts for a cell to be included in the output matrix
    'min_features' : the minimum number of features for a cell to be included in the output matrix
    'min_cells' : the minimum number of cells for a feature to be included in the output matrix
    'l_features' : a list of features (values in the first column of 'features.tsv.gz') to include. All other features will be excluded from the output matrix. (default: None) If None is given, include all features in the output matrix.
    'l_cells' : a list of cells (values in the first column of 'barcodes.tsv.gz') to include. All other cells will be excluded from the output matrix. (default: None) If None is given, include all cells in the output matrix.
    'int_num_threads' : when 'int_num_threads' is 1, does not use the multiprocessing  module for parallel processing
    'function_for_adjusting_thresholds' : a function for adjusting thresholds based on the summarized metrics. Useful when the exact threshold for removing empty droplets are variable across the samples. the function should receive arguments and return values in the following structure:
                                        min_counts_new, min_features_new, min_cells_new = function_for_adjusting_thresholds( path_folder_mtx_10x_output, min_counts, min_features, min_cells )
    """

    """ handle inputs """
    if path_folder_mtx_10x_input[-1] != "/":
        path_folder_mtx_10x_input += "/"
    if path_folder_mtx_10x_output[-1] != "/":
        path_folder_mtx_10x_output += "/"
    if ((min_counts is not None) or (min_features is not None)) and (
        min_cells is not None
    ):  # check whether thresholds for both cells and features were given (thresdholds for either cells or features can be given at a time)
        if verbose:
            logger.info(
                "[MTX_10X_Filter] (error) no threshold is given or more thresholds for both cells and features are given. (Thresdholds for either cells or features can be given at a time.)"
            )
        return -1
    # create an output folder
    filesystem_operations("mkdir", path_folder_mtx_10x_output, exist_ok=True)

    # define input file directories
    path_file_input_bc = f"{path_folder_mtx_10x_input}barcodes.tsv.gz"
    path_file_input_feature = f"{path_folder_mtx_10x_input}features.tsv.gz"
    path_file_input_mtx = f"{path_folder_mtx_10x_input}matrix.mtx.gz"

    str_datatype = _MTX_Detect_data_type(
        path_file_input_mtx
    )  # retrieve the data type of the matrix in string format

    # check whether all required files are present
    if sum(
        list(
            not filesystem_operations("exists", path_folder)
            for path_folder in [
                path_file_input_bc,
                path_file_input_feature,
                path_file_input_mtx,
            ]
        )
    ):
        if verbose:
            logger.info(f"required file(s) is not present at {path_folder_mtx_10x}")

    """ read barcode and feature information """
    df_bc = pd.read_csv(path_file_input_bc, sep="\t", header=None)
    df_bc.columns = ["barcode"]
    df_feature = pd.read_csv(path_file_input_feature, sep="\t", header=None)
    df_feature.columns = ["id_feature", "feature", "feature_type"]

    """ split input mtx file into multiple files """
    l_path_file_mtx_10x = MTX_10X_Split(
        path_folder_mtx_10x_input,
        int_max_num_entries_for_chunk=int_max_num_entries_for_chunk,
        flag_split_mtx=flag_split_mtx,
    )

    """ summarizes counts """
    dict_data = MTX_10X_Summarize_Counts(
        path_folder_mtx_10x_input,
        verbose=verbose,
        int_num_threads=int_num_threads,
        flag_split_mtx=flag_split_mtx,
        int_max_num_entries_for_chunk=int_max_num_entries_for_chunk,
    )
    (
        dict_id_column_to_count,
        dict_id_column_to_n_features,
        dict_id_row_to_count,
        dict_id_row_to_n_cells,
        dict_id_row_to_log_transformed_count,
    ) = (
        dict_data["dict_id_column_to_count"],
        dict_data["dict_id_column_to_n_features"],
        dict_data["dict_id_row_to_count"],
        dict_data["dict_id_row_to_n_cells"],
        dict_data["dict_id_row_to_log_transformed_count"],
    )  # parse 'dict_data'

    """ adjust thresholds based on the summarized metrices (if a function has been given) """
    if function_for_adjusting_thresholds is not None:
        min_counts, min_features, min_cells = function_for_adjusting_thresholds(
            path_folder_mtx_10x_input, min_counts, min_features, min_cells
        )

    """ filter row or column that do not satisfy the given thresholds """
    if min_counts is not None:
        dict_id_column_to_count = dict(
            (k, dict_id_column_to_count[k])
            for k in dict_id_column_to_count
            if dict_id_column_to_count[k] >= min_counts
        )
    if min_features is not None:
        dict_id_column_to_n_features = dict(
            (k, dict_id_column_to_n_features[k])
            for k in dict_id_column_to_n_features
            if dict_id_column_to_n_features[k] >= min_features
        )
    if min_cells is not None:
        dict_id_row_to_n_cells = dict(
            (k, dict_id_row_to_n_cells[k])
            for k in dict_id_row_to_n_cells
            if dict_id_row_to_n_cells[k] >= min_cells
        )

    """ retrieve id_row and id_column that satisfy the given thresholds """
    set_id_column = set(dict_id_column_to_count).intersection(
        set(dict_id_column_to_n_features)
    )
    set_id_row = set(dict_id_row_to_n_cells)

    """ exclude cells and features not present in the input lists (if the lists were given)  """
    if l_cells is not None:
        dict_barcode_to_id_column = dict(
            (barcode, id_column)
            for id_column, barcode in enumerate(df_bc.barcode.values)
        )
        set_id_column = set_id_column.intersection(
            set(
                dict_barcode_to_id_column[barcode]
                for barcode in set(l_cells)
                if barcode in dict_barcode_to_id_column
            )
        )
        del dict_barcode_to_id_column
    if l_features is not None:
        dict_id_feature_to_id_row = dict(
            (id_feature, id_row)
            for id_row, id_feature in enumerate(df_feature.id_feature.values)
        )
        set_id_row = set_id_row.intersection(
            set(
                dict_id_feature_to_id_row[id_feature]
                for id_feature in set(l_features)
                if id_feature in dict_id_feature_to_id_row
            )
        )
        del dict_id_feature_to_id_row

    """ report the number of cells or features that will be filtered out """
    if verbose:
        int_n_bc_filtered = len(df_bc) - len(set_id_column)
        if int_n_bc_filtered > 0:
            logger.info(
                f"{int_n_bc_filtered}/{len( df_bc )} barcodes will be filtered out"
            )
        int_n_feature_filtered = len(df_feature) - len(set_id_row)
        if int_n_feature_filtered > 0:
            logger.info(
                f"{int_n_feature_filtered}/{len( df_feature )} features will be filtered out"
            )

    """ retrieve a mapping between previous id_column to current id_column """
    global dict_id_column_previous_to_id_column_current, dict_id_row_previous_to_id_row_current  # use global variables for multiprocessing
    df_bc = df_bc.loc[list(set_id_column)]
    df_bc.index.name = "id_column_previous"
    df_bc.reset_index(drop=False, inplace=True)
    df_bc["id_column_current"] = np.arange(len(df_bc))
    dict_id_column_previous_to_id_column_current = df_bc.set_index(
        "id_column_previous"
    ).id_column_current.to_dict()
    bk.PICKLE_Write(
        f"{path_folder_mtx_10x_output}dict_id_column_previous_to_id_column_current.pickle",
        dict_id_column_previous_to_id_column_current,
    )  # save id_feature to index_feature mapping
    """ retrieve a mapping between previous id_row to current id_row """
    df_feature = df_feature.loc[list(set_id_row)]
    df_feature.index.name = "id_row_previous"
    df_feature.reset_index(drop=False, inplace=True)
    df_feature["id_row_current"] = np.arange(len(df_feature))
    dict_id_row_previous_to_id_row_current = df_feature.set_index(
        "id_row_previous"
    ).id_row_current.to_dict()
    bk.PICKLE_Write(
        f"{path_folder_mtx_10x_output}dict_id_row_previous_to_id_row_current.pickle",
        dict_id_row_previous_to_id_row_current,
    )  # save id_feature to index_feature mapping

    """ save barcode file """
    df_bc.to_csv(
        f"{path_folder_mtx_10x_output}barcodes.tsv.gz",
        columns=["barcode"],
        sep="\t",
        index=False,
        header=False,
    )
    del df_bc

    """ save feature file """
    df_feature[["id_feature", "feature", "feature_type"]].to_csv(
        f"{path_folder_mtx_10x_output}features.tsv.gz",
        sep="\t",
        index=False,
        header=False,
    )  # save as a file
    del df_feature

    """ write a filtered matrix.mtx.gz for each split mtx file using multiple processes and retrieve the total number of entries written by each process """
    # compose inputs for multiprocessing
    df_input = pd.DataFrame(
        {
            "path_file_mtx_10x": l_path_file_mtx_10x,
            "index_mtx_10x": np.arange(len(l_path_file_mtx_10x)),
        }
    )

    l_int_n_entries = bk.Multiprocessing(
        df_input,
        __MTX_10X_Filter__filter_mtx_10x__,
        int_num_threads,
        global_arguments=[path_folder_mtx_10x_output, str_datatype],
    )
    # retrieve the total number of entries
    int_total_n_entries = sum(l_int_n_entries)

    """ combine parts and add the MTX file header to compose a combined mtx file """
    df_file = bk.GLOB_Retrive_Strings_in_Wildcards(
        f"{path_folder_mtx_10x_output}matrix.mtx.gz.*.gz"
    )
    df_file.wildcard_0 = df_file.wildcard_0.astype(int)
    df_file.sort_values("wildcard_0", inplace=True)

    # write header
    path_file_header = f"{path_folder_mtx_10x_output}matrix.mtx.header.txt.gz"
    with gzip.open(path_file_header, "wb") as newfile:
        newfile.write(
            f"%%MatrixMarket matrix coordinate {str_datatype} general\n%\n{len( dict_id_row_previous_to_id_row_current )} {len( dict_id_column_previous_to_id_column_current )} {int_total_n_entries}\n".encode()
        )
    bk.OS_Run(
        ["cat", path_file_header] + list(df_file.path.values),
        path_file_stdout=f"{path_folder_mtx_10x_output}matrix.mtx.gz",
        stdout_binary=True,
        return_output=False,
    )  # combine the output mtx files in the order # does not delete temporary files if 'flag_split_mtx' is True

    # write a flag indicating that the current output directory contains split mtx files
    with open(f"{path_folder_mtx_10x_output}matrix.mtx.gz.split.flag", "w") as file:
        file.write("completed")


def MTX_10X_Identify_Highly_Variable_Features(
    path_folder_mtx_10x_input,
    int_target_sum=10000,
    verbose=False,
    int_num_threads=15,
    flag_split_mtx=True,
    int_max_num_entries_for_chunk=10000000,
):
    """# 2022-03-16 17:18:44
    calculate variance from log-transformed normalized counts using 'MTX_10X_Calculate_Average_Log10_Transformed_Normalized_Expr' and rank features based on how each feature is variable compared to other features with similar means.
    Specifically, polynomial of degree 2 will be fitted to variance-mean relationship graph in order to captures the relationship between variance and mean.

    'name_col_for_mean', 'name_col_for_variance' : name of columns of 'df_summary' returned by 'MTX_10X_Calculate_Average_Log10_Transformed_Normalized_Expr' that will be used to infer highly variable features. By defaults, mean and variance of log-transformed normalized counts will be used.
    """

    # calculate variance and means and load the result
    df_summary = MTX_10X_Calculate_Average_Log10_Transformed_Normalized_Expr(
        path_folder_mtx_10x_input,
        int_target_sum=int_target_sum,
        int_num_threads=int_num_threads,
        verbose=verbose,
        flag_split_mtx=flag_split_mtx,
    )

    # calculate scores for identifying highly variable features for the selected set of count data types: [ 'log_transformed_normalized_count', 'log_transformed_count' ]
    for name_type in ["log_transformed_normalized_count", "log_transformed_count"]:
        name_col_for_mean, name_col_for_variance = (
            f"mean_{name_type}",
            f"variance_of_{name_type}",
        )
        # retrieve the relationship between mean and variance
        arr_mean = df_summary[name_col_for_mean].values
        arr_var = df_summary[name_col_for_variance].values
        mean_var_relationship_fit = np.polynomial.polynomial.Polynomial.fit(
            arr_mean, arr_var, 2
        )

        # calculate the deviation from the estimated variance from the mean
        arr_ratio_of_variance_to_expected_variance_from_mean = np.zeros(len(df_summary))
        arr_diff_of_variance_to_expected_variance_from_mean = np.zeros(len(df_summary))
        for i in range(len(df_summary)):  # iterate list of means of the features
            var, mean = arr_var[i], arr_mean[i]  # retrieve var and mean
            var_expected = mean_var_relationship_fit(
                mean
            )  # calculate expected variance from the mean
            if (
                var_expected == 0
            ):  # handle the case when the current expected variance is zero
                arr_ratio_of_variance_to_expected_variance_from_mean[i] = 1
                arr_diff_of_variance_to_expected_variance_from_mean[i] = 0
            else:
                arr_ratio_of_variance_to_expected_variance_from_mean[i] = (
                    var / var_expected
                )
                arr_diff_of_variance_to_expected_variance_from_mean[i] = (
                    var - var_expected
                )

        df_summary[
            f"float_ratio_of_variance_to_expected_variance_from_mean_from_{name_type}"
        ] = arr_ratio_of_variance_to_expected_variance_from_mean
        df_summary[
            f"float_diff_of_variance_to_expected_variance_from_mean_{name_type}"
        ] = arr_diff_of_variance_to_expected_variance_from_mean

        # calculate the product of the ratio and difference of variance to expected variance for scoring and sorting highly variable features
        df_summary[f"float_score_highly_variable_feature_from_{name_type}"] = (
            df_summary[
                f"float_ratio_of_variance_to_expected_variance_from_mean_from_{name_type}"
            ]
            * df_summary[
                f"float_diff_of_variance_to_expected_variance_from_mean_{name_type}"
            ]
        )

    df_summary["float_score_highly_variable_feature"] = list(
        np.prod(arr_val) if sum(np.sign(arr_val) < 0) == 0 else 0
        for arr_val in df_summary[
            [
                "float_score_highly_variable_feature_from_log_transformed_normalized_count",
                "float_score_highly_variable_feature_from_log_transformed_count",
            ]
        ].values
    )
    return df_summary


""" functions that were added later """


def is_binary_stream(f):
    """# 2022-05-01 01:57:10
    check whether a given stream is a binary stream"""
    if hasattr(f, "mode"):  # if given stream is file
        return "b" in f.mode
    else:
        return isinstance(f, (io.RawIOBase, io.BufferedIOBase))


def __Get_path_essential_files__(path_folder_mtx_10x_input):
    """# 2022-04-30 16:28:15
    get paths of essential files for the given matrix folder ('path_folder_mtx_10x_input', currently only supports 10X-GEX-formated output matrix)
    """
    # define input file paths
    path_file_input_bc = f"{path_folder_mtx_10x_input}barcodes.tsv.gz"
    path_file_input_feature = f"{path_folder_mtx_10x_input}features.tsv.gz"
    path_file_input_mtx = f"{path_folder_mtx_10x_input}matrix.mtx.gz"
    # check whether input files exist
    for path_file in [path_file_input_bc, path_file_input_feature, path_file_input_mtx]:
        if not os.path.exists(path_file):
            raise OSError(f"{path_file} does not exist")
    return path_file_input_bc, path_file_input_feature, path_file_input_mtx


def Merge_Sort_Files(file_output, *l_iterator_decorated_file_input):
    """# 2022-05-01 02:23:09
    Merge sort input files (should be sorted) without loading the complete contents on memory.
    'path_file_output' : output file handle/stream
    'l_iterator_decorated_file_input' : a list of iterators based on input file handles (or streams). each iterator should yield the following tuple: (key_for_sorting, content_that_will_be_written_in_the_output_file). This function does not check whether the datatype of the 'content_that_will_be_written_in_the_output_file' matches that of 'path_file_output'
    """
    # handle invalid case
    if len(l_iterator_decorated_file_input) == 0:
        return -1
    # perform merge sorting
    for r in heapq.merge(*l_iterator_decorated_file_input):
        file_output.write(
            r[1]
        )  # assumes the 'iterator_decorated_file_input' returns appropriate datatype (either bytes or string) for the output file


def __Merge_Sort_MTX_10X__(
    path_file_output,
    *l_path_file_input,
    flag_ramtx_sorted_by_id_feature=True,
    flag_delete_input_file_upon_completion=False,
):
    """# 2022-05-01 02:25:07
    merge sort mtx files
    'path_file_output' and 'l_path_file_input'  : either file path or file handles

    'flag_ramtx_sorted_by_id_feature' : if True, sort by 'id_feature'. if False, sort by 'id_cell'
    """
    # process arguments for input files
    if isinstance(l_path_file_input[0], str):  # if paths are given as input files
        flag_input_binary = (
            l_path_file_input[0].rsplit(".", 1)[1].lower() == "gz"
        )  # automatically detect gzipped input file # determined gzipped status by only looking at the first file
        l_file_input = list(
            gzip.open(path_file, "rb") if flag_input_binary else open(path_file, "r")
            for path_file in l_path_file_input
        )
    else:
        flag_input_binary = is_binary_stream(l_file_input[0])  # detect binary stream
        l_file_input = l_path_file_input
    # process argument for output file
    if isinstance(path_file_output, str):  # if path was given as an output file
        flag_output_is_file = True
        flag_output_binary = (
            path_file_output.rsplit(".", 1)[1].lower() == "gz"
        )  # automatically detect gzipped input file # determined gzipped status by only looking at the first file
        file_output = (
            gzip.open(path_file_output, "wb")
            if flag_output_binary
            else open(path_file_output, "w")
        )
    else:
        flag_output_is_file = False
        flag_output_binary = is_binary_stream(path_file_output)  # detect binary stream
        file_output = path_file_output

    # define a function for decorating mtx record
    def __decorate_mtx_file__(file):
        while True:
            line = file.readline()
            if len(line) == 0:
                break
            """ parse a mtx record """
            line_decoded = line.decode() if flag_input_binary else line
            index_row, index_column, float_value = (
                (line_decoded).strip().split()
            )  # parse a record of a matrix-market format file
            index_row, index_column, float_value = (
                int(index_row),
                int(index_column),
                float(float_value),
            )  # 0-based coordinates
            yield index_row if flag_ramtx_sorted_by_id_feature else index_column, (
                (line if flag_input_binary else line.encode())
                if flag_output_binary
                else line_decoded
            )

    Merge_Sort_Files(
        file_output, *list(__decorate_mtx_file__(file) for file in l_file_input)
    )  # perform merge sorting

    # if the output file is stream, does not close the stream # only close the output if the output file was an actual file
    if flag_output_is_file:
        file_output.close()

    """ delete input files once merge sort is completed if 'flag_delete_input_file_upon_completion' is True """
    if flag_delete_input_file_upon_completion and isinstance(
        l_path_file_input[0], str
    ):  # if paths are given as input files
        for path_file in l_path_file_input:
            os.remove(path_file)


def __Merge_Sort_and_Index_MTX_10X__(
    path_file_output,
    *l_path_file_input,
    flag_ramtx_sorted_by_id_feature=True,
    flag_delete_input_file_upon_completion=False,
):
    """# 2022-05-01 02:25:07
    merge sort mtx files into a single mtx uncompressed file and index entries in the combined mtx file while writing the file
    'path_file_output' : should be a file path, file handle (or stream) for non-binary (text) output
    'l_path_file_input'

    'flag_ramtx_sorted_by_id_feature' : if True, sort by 'id_feature'. if False, sort by 'id_cell'
    """
    # process arguments for input files
    if isinstance(l_path_file_input[0], str):  # if paths are given as input files
        flag_input_binary = (
            l_path_file_input[0].rsplit(".", 1)[1].lower() == "gz"
        )  # automatically detect gzipped input file # determined gzipped status by only looking at the first file
        l_file_input = list(
            gzip.open(path_file, "rb") if flag_input_binary else open(path_file, "r")
            for path_file in l_path_file_input
        )
    else:
        flag_input_binary = is_binary_stream(l_file_input[0])  # detect binary stream
        l_file_input = l_path_file_input
    # process argument for output file
    if isinstance(path_file_output, str):  # if path was given as an output file
        flag_output_is_file = True
        flag_output_binary = (
            path_file_output.rsplit(".", 1)[1].lower() == "gz"
        )  # automatically detect gzipped input file # determined gzipped status by only looking at the first file
        file_output = (
            gzip.open(path_file_output, "wb")
            if flag_output_binary
            else open(path_file_output, "w")
        )
    else:
        flag_output_is_file = False
        flag_output_binary = is_binary_stream(path_file_output)  # detect binary stream
        file_output = path_file_output

    if flag_output_binary:  # the output file should be non-binary stream/file
        raise OSError("the output file should be non-binary stream/file")

    # define and open index output file
    path_file_index_output = f"{path_file_output}.idx.tsv.gz"
    file_index_output = gzip.open(path_file_index_output, "wb")
    file_index_output.write(
        ("\t".join(["index_entry", "int_pos_start", "int_pos_end"]) + "\n").encode()
    )  # write the header of the index file

    # define a function for decorating mtx record
    def __decorate_mtx_file__(file):
        while True:
            line = file.readline()
            if len(line) == 0:
                break
            """ parse a mtx record """
            line_decoded = line.decode() if flag_input_binary else line
            index_row, index_column, float_value = (
                (line_decoded).strip().split()
            )  # parse a record of a matrix-market format file
            index_row, index_column, float_value = (
                int(index_row),
                int(index_column),
                float(float_value),
            )  # 0-based coordinates
            yield index_row if flag_ramtx_sorted_by_id_feature else index_column, (
                (line if flag_input_binary else line.encode())
                if flag_output_binary
                else line_decoded
            )

    # perform merge sorting
    index_entry_currently_being_written = -1
    int_num_character_written_for_index_entry_currently_being_written = 0
    int_total_num_character_written = 0
    for r in heapq.merge(*list(__decorate_mtx_file__(file) for file in l_file_input)):
        if (
            index_entry_currently_being_written != r[0]
        ):  # if current index_entry is different from the previous one, which mark the change of sorted blocks (a block has the same id_entry), record the data for the previous block and initialze data for the next block
            if (
                index_entry_currently_being_written > 0
            ):  # check whether 'index_entry_currently_being_written' is valid (ignore 'dummy' or default value that was used for initialization)
                file_index_output.write(
                    (
                        "\t".join(
                            map(
                                str,
                                [
                                    index_entry_currently_being_written,
                                    int_total_num_character_written,
                                    int_total_num_character_written
                                    + int_num_character_written_for_index_entry_currently_being_written,
                                ],
                            )
                        )
                        + "\n"
                    ).encode()
                )  # write information required for indexing
            int_total_num_character_written += int_num_character_written_for_index_entry_currently_being_written  # update 'int_total_num_character_written'
            # initialize data for index of the next 'index_entry'
            index_entry_currently_being_written = r[0]  # update current index_entry
            int_num_character_written_for_index_entry_currently_being_written = 0  # reset the count of characters (which is equal to the number of bytes for any mtx records, because they only contains numeric characters)
        int_num_character_written_for_index_entry_currently_being_written += file_output.write(
            r[1]
        )  # assumes the 'iterator_decorated_file_input' returns appropriate datatype (either bytes or string) for the output file # count the number of characters written for the current index_row

    # write the record for the last block
    file_index_output.write(
        (
            "\t".join(
                map(
                    str,
                    [
                        index_entry_currently_being_written,
                        int_total_num_character_written,
                        int_total_num_character_written
                        + int_num_character_written_for_index_entry_currently_being_written,
                    ],
                )
            )
            + "\n"
        ).encode()
    )  # write information required for indexing
    # close index file
    file_index_output.close()
    # if the output file is stream, does not close the stream # only close the output if the output file was an actual file
    if flag_output_is_file:
        file_output.close()

    """ delete input files once merge sort is completed if 'flag_delete_input_file_upon_completion' is True """
    if flag_delete_input_file_upon_completion and isinstance(
        l_path_file_input[0], str
    ):  # if paths are given as input files
        for path_file in l_path_file_input:
            os.remove(path_file)


""" methods for handling 10X matrix objects """


def Convert_df_count_to_MTX_10X(
    path_file_df_count: str,
    path_folder_mtx_10x_output: str,
    path_folder_mtx_10x_filtered_output: str,
    chunksize: int = 1000000,
    int_min_count_features_for_filtering_barcodes: int = 50,
    flag_output_dtype_is_integer: bool = True,
):
    """# 2023-01-06 23:46:20
    convert df_count (ouro output) to 10X MTX (matrix market) format in a memory-efficient manner.

    path_file_df_count : str, # file path to 'df_count'
    path_folder_mtx_10x_output : str, # a folder containing 10x output matrix (unfiltered)
    path_folder_mtx_10x_filtered_output : str, # a folder containing 10x output matrix (filtered)
    chunksize : int = 500000,
    int_min_count_features_for_filtering_barcodes : int = 50, # the minimum number of features in a barcode to be included in the filtered output
    flag_output_dtype_is_integer : bool = True, # a boolean flag indicating the output dtype is integer dtype. Set this flag to False if the output dtype is float
    """
    # create output folders
    os.makedirs(path_folder_mtx_10x_output, exist_ok=True)
    os.makedirs(path_folder_mtx_10x_filtered_output, exist_ok=True)

    # retrieve unique feature/barcode information from df_count without loading entire data in the memory
    bk.DF_Deduplicate_without_loading_in_memory(
        path_file_df_count,
        f"{path_folder_mtx_10x_output}_features.tsv.gz",
        l_col_for_identifying_duplicates=["feature", "id_feature"],
        str_delimiter="\t",
    )
    res = bk.DF_Deduplicate_without_loading_in_memory(
        path_file_df_count,
        f"{path_folder_mtx_10x_output}_barcodes.tsv.gz",
        l_col_for_identifying_duplicates="barcode",
        str_delimiter="\t",
    )  # collect the number of records
    int_num_lines = res["int_num_lines"]
    s_num_records_for_each_barcode = pd.Series(
        res["dict_t_val_count"]
    )  # retrieve the number of records for each barcode
    del res
    s_num_records_for_each_barcode = s_num_records_for_each_barcode[
        s_num_records_for_each_barcode >= int_min_count_features_for_filtering_barcodes
    ]  # filter barcodes using the given setting
    df_barcode_filtered = pd.DataFrame(
        {"barcode": s_num_records_for_each_barcode.index.values}
    )  # compose a dataframe containing filtered barcodes

    # read features and barcode information
    df_barcode = pd.read_csv(
        f"{path_folder_mtx_10x_output}_barcodes.tsv.gz", sep="\t", usecols=["barcode"]
    )
    df_feature = pd.read_csv(
        f"{path_folder_mtx_10x_output}_features.tsv.gz",
        sep="\t",
        usecols=["feature", "id_feature"],
    )
    df_feature = df_feature.loc[:, ["id_feature", "feature"]]
    df_feature["10X_type"] = "Gene Expression"
    # save feature/cell metadata
    df_barcode.to_csv(
        f"{path_folder_mtx_10x_output}barcodes.tsv.gz",
        sep="\t",
        index=False,
        header=False,
    )
    df_barcode_filtered.to_csv(
        f"{path_folder_mtx_10x_filtered_output}barcodes.tsv.gz",
        sep="\t",
        index=False,
        header=False,
    )
    df_feature.to_csv(
        f"{path_folder_mtx_10x_output}features.tsv.gz",
        sep="\t",
        index=False,
        header=False,
    )
    df_feature.to_csv(
        f"{path_folder_mtx_10x_filtered_output}features.tsv.gz",
        sep="\t",
        index=False,
        header=False,
    )

    # retrieve barcode/feature to integer representation of barcode/feature mapping
    dict_to_int_barcode = dict(
        (e, i + 1) for i, e in enumerate(df_barcode.iloc[:, 0].values)
    )
    dict_to_int_barcode_filtered = dict(
        (e, i + 1) for i, e in enumerate(df_barcode_filtered.iloc[:, 0].values)
    )
    dict_to_int_feature = dict(
        (e, i + 1) for i, e in enumerate(df_feature.iloc[:, 0].values)
    )

    (
        int_num_features,
        int_num_barcodes,
        int_num_barcodes_filtered,
        int_num_records,
        int_num_records_filtered,
    ) = (
        len(df_feature),
        len(df_barcode),
        len(df_barcode_filtered),
        int_num_lines,
        s_num_records_for_each_barcode.sum(),
    )  # retrieve metadata of the output matrix
    del (
        df_feature,
        df_barcode,
        df_barcode_filtered,
        s_num_records_for_each_barcode,
    )  # delete objects

    # write mtx file
    str_dtype_mtx_header = (
        "integer" if flag_output_dtype_is_integer else "real"
    )  # retrieve dtype name for the matrix header
    with gzip.open(f"{path_folder_mtx_10x_output}matrix.mtx.gz", "wb") as newfile:
        newfile.write(
            f"""%%MatrixMarket matrix coordinate {str_dtype_mtx_header} general\n%\n{int_num_features} {int_num_barcodes} {int_num_records}\n""".encode()
        )
        with gzip.open(
            f"{path_folder_mtx_10x_filtered_output}matrix.mtx.gz", "wb"
        ) as newfile_filtered:
            newfile_filtered.write(
                f"""%%MatrixMarket matrix coordinate {str_dtype_mtx_header} general\n%\n{int_num_features} {int_num_barcodes_filtered} {int_num_records_filtered}\n""".encode()
            )

            # iterate through each chunk
            for df_chunk in pd.read_csv(
                path_file_df_count,
                iterator=True,
                header=0,
                chunksize=chunksize,
                sep="\t",
                usecols=["id_feature", "barcode", "read_count"],
            ):
                df_chunk = df_chunk[
                    ["id_feature", "barcode", "read_count"]
                ]  # reorder columns
                mask_filtered = np.array(
                    list(
                        e in dict_to_int_barcode_filtered
                        for e in df_chunk.barcode.values
                    ),
                    dtype=bool,
                )  # retrieve a mask for filtering records
                df_chunk["id_feature"] = df_chunk.id_feature.apply(
                    MAP.Map(dict_to_int_feature).a2b
                )
                df_chunk_filtered = df_chunk[mask_filtered]  # filter records
                df_chunk["barcode"] = df_chunk.barcode.apply(
                    MAP.Map(dict_to_int_barcode).a2b
                )
                df_chunk_filtered["barcode"] = df_chunk_filtered.barcode.apply(
                    MAP.Map(dict_to_int_barcode_filtered).a2b
                )
                df_chunk.to_csv(newfile, sep=" ", header=None, index=False)
                df_chunk_filtered.to_csv(
                    newfile_filtered, sep=" ", header=None, index=False
                )  # write filtered records
                del mask_filtered, df_chunk, df_chunk_filtered

    # delete temporary files
    os.remove(f"{path_folder_mtx_10x_output}_barcodes.tsv.gz")
    os.remove(f"{path_folder_mtx_10x_output}_features.tsv.gz")


def MTX_Convert_10x_MEX_to_10x_HDF5_Format(
    path_folder_matrix_input_mex_format: str,
    path_file_matrix_output_hdf5_format: str,
    name_genome: str = "unknown",
):
    """# 2023-09-15 01:36:49
    path_folder_matrix_input_mex_format : str # the path of the input 10x MEX matrix folder
    path_file_matrix_output_hdf5_format : str # the path of the output 10x HDF5 matrix file
    name_genome : str = 'unknown' # the name of the genome
    """
    """ import libaries """
    import h5py

    """ read 10x MEX format """
    # read mtx file as a tabular format
    df_mtx = pd.read_csv(
        f"{path_folder_matrix_input_mex_format}matrix.mtx.gz", sep=" ", comment="%"
    )
    df_mtx.columns = ["id_row", "id_column", "read_count"]
    df_mtx.sort_values("id_column", inplace=True)  # sort using id_cell
    # read barcodes
    arr_bc = pd.read_csv(
        f"{path_folder_matrix_input_mex_format}barcodes.tsv.gz", sep="\t", header=None
    ).values.ravel()
    # read feature tables
    df_feature = pd.read_csv(
        f"{path_folder_matrix_input_mex_format}features.tsv.gz", sep="\t", header=None
    )
    df_feature.columns = ["id_feature", "feature", "feature_type"]

    """ write hdf5 file """
    newfile = h5py.File(path_file_matrix_output_hdf5_format, "w")  # open new HDF5 file

    def _write_string_array(handle, name_array: str, arr_str: List[str]):
        """# 2023-09-14 21:41:14
        write a string array to a HDF5 object
        """
        handle.create_dataset(
            name_array,
            (len(arr_str),),
            dtype="S" + str(np.max(list(len(e) for e in arr_str))),
            data=list(e.encode("ascii", "ignore") for e in arr_str),
        )  # writing string dtype array

    # create matrix group
    mtx = newfile.create_group("matrix")

    # write barcodes
    _write_string_array(mtx, "barcodes", arr_bc)

    # # write id/names

    # write data
    arr = df_mtx.read_count.values
    flag_dtype_is_integer = np.issubdtype(arr.dtype, np.integer)  # check integer dtype
    mtx.create_dataset("data", (len(arr),), "i8" if flag_dtype_is_integer else "f", arr)

    # write indices
    arr = df_mtx.id_row.values - 1  # 1 -> 0-based coordinates
    mtx.create_dataset("indices", (len(arr),), "i8", arr)

    # write shape
    mtx.create_dataset("shape", (2,), "i8", [len(df_feature), len(arr_bc)])

    # write indptr
    arr = df_mtx.id_column.values
    arr = arr - 1  # 1>0-based coordinate
    int_num_bc = len(arr_bc)  # retrieve the number of barcodes
    int_num_records = len(arr)  # retrieve the number of records
    arr_indptr = np.zeros(int_num_bc + 1, dtype="i8")  # initialize 'arr_indptr'
    arr_indptr[-1] = int_num_records  # last entry should be the number of the records
    id_col_current = arr[0]  # initialize 'id_col_current'
    for i, id_col in enumerate(arr):
        if id_col_current != id_col:
            if (
                id_col_current + 1 < id_col
            ):  # if there are some skipped columns ('barcodes' with zero number of records)
                for id_col_with_no_records in range(id_col_current + 1, id_col):
                    arr_indptr[id_col_with_no_records] = (
                        i  # add 'indptr' for the 'barcodes' with zero number of records
                    )
            id_col_current = id_col  # update 'id_col_current'
            arr_indptr[id_col] = i
    if id_col_current + 1 < int_num_bc:
        for id_col_with_no_records in range(id_col_current + 1, int_num_bc):
            arr_indptr[id_col_with_no_records] = (
                int_num_records  # add 'indptr' for the 'barcodes' with zero number of records
            )
    mtx.create_dataset("indptr", (len(arr_indptr),), "i8", arr_indptr)

    # create matrix group
    ft = mtx.create_group("features")

    # write features/id, features/name, features/feature_type
    _write_string_array(ft, "id", df_feature.id_feature.values)
    _write_string_array(ft, "name", df_feature.feature.values)
    _write_string_array(ft, "feature_type", df_feature.feature_type.values)
    _write_string_array(
        ft, "genome", [name_genome] * len(df_feature)
    )  # add genome data type (required for scanpy)

    # close the file
    newfile.close()


def MTX_Convert_10x_HDF5_to_10x_MEX_Format(
    path_file_matrix_input_hdf5_format: str,
    path_folder_matrix_output_mex_format: str,
):
    """# 2023-10-04 13:31:40
    path_file_matrix_input_hdf5_format: str, # the path of the input 10x HDF5 matrix file
    path_folder_matrix_output_mex_format: str, # the path of the output 10x MEX matrix folder
    """
    """ import libaries """
    import h5py
    import scipy.io

    """ read hdf5 file """
    newfile = h5py.File(path_file_matrix_input_hdf5_format, "r")  # open new HDF5 file
    mtx = newfile["matrix"]  # retrieve the group

    """ read barcodes """
    arr_bc = mtx["barcodes"][:]  # retrieve col (cell) boundary positions
    arr_bc = np.array(
        list(e.decode() for e in arr_bc), dtype=object
    )  # change dtype of the barcode

    """ read features """
    ft = mtx["features"]  # retrieve the group
    arr_id_ft = np.array(
        list(e.decode() for e in ft["id"][:]), dtype=object
    )  # change dtype of the barcode
    arr_id_name = np.array(
        list(e.decode() for e in ft["name"][:]), dtype=object
    )  # change dtype of the barcode
    arr_id_feature_type = np.array(
        list(e.decode() for e in ft["feature_type"][:]), dtype=object
    )  # change dtype of the barcode
    arr_genome = np.array(
        list(e.decode() for e in ft["genome"][:]), dtype=object
    )  # change dtype of the barcode
    # compose feature dataframe
    df_feature = pd.DataFrame(
        {
            "id_feature": arr_id_ft,
            "feature": arr_id_name,
            "feature_type": arr_id_feature_type,
            "genome": arr_genome,
        }
    )

    """ read count values """
    arr_data = mtx["data"][:]
    arr_row_indices = mtx["indices"][:]  # retrieve row (gene) indices
    arr_col_index_boundary = mtx["indptr"][:]  # retrieve col (cell) boundary positions
    # compose arr_col_indices
    arr_col_indices = np.zeros_like(arr_row_indices)  # initialize 'arr_col_indices'
    for idx_bc in range(len(arr_bc)):  # for each barcode index (0-based coordinates)
        arr_col_indices[
            arr_col_index_boundary[idx_bc] : arr_col_index_boundary[idx_bc + 1]
        ] = idx_bc
    # compose 'df_mtx'
    df_mtx = pd.DataFrame(
        {
            "id_row": arr_row_indices,  # 0 based coordinates
            "id_column": arr_col_indices,  # 0 based coordinates
            "read_count": arr_data,
        }
    )

    """ write the output MEX files """
    # create the output directory
    os.makedirs(path_folder_matrix_output_mex_format, exist_ok=True)

    # save barcodes
    pd.DataFrame(arr_bc).to_csv(
        f"{path_folder_matrix_output_mex_format}barcodes.tsv.gz",
        sep="\t",
        index=False,
        header=False,
    )

    # save features
    # compose a feature dataframe
    df_feature[["id_feature", "feature", "feature_type"]].to_csv(
        f"{path_folder_matrix_output_mex_format}features.tsv.gz",
        sep="\t",
        index=False,
        header=False,
    )  # save as a file
    # retrieve list of features

    # save count matrix as a gzipped matrix market format
    row, col, data = df_mtx[["id_row", "id_column", "read_count"]].values.T
    sm = scipy.sparse.coo_matrix(
        (data, (row, col)), shape=(len(df_feature), len(arr_bc))
    )
    scipy.io.mmwrite(f"{path_folder_matrix_output_mex_format}matrix", sm)

    # remove previous output file to overwrite the file
    path_file_mtx_output = f"{path_folder_matrix_output_mex_format}matrix.mtx.gz"
    if os.path.exists(path_file_mtx_output):
        os.remove(path_file_mtx_output)
    bk.OS_Run(
        ["gzip", f"{path_folder_matrix_output_mex_format}matrix.mtx"]
    )  # gzip the mtx file
