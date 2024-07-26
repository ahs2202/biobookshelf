"""
classes and functions for sharing data across multiple forked processes
"""

from multiprocessing.managers import BaseManager
import numpy as np
from typing import Union, List, Dict
import zarr


class ManagerFileSystem(BaseManager):
    pass


class FileSystemOperator:
    """# 2023-09-24 14:50:47
    A class intended for performing asynchronous file system operations in a separate, managed process. By using multiple managers, concurrent, asynchronous operations can be performed in multiple processes. These managers can be used multiple times.

    dict_kwargs_s3 : dict = dict( ) # s3 credentials to use
    """

    # constructor
    def __init__(self, dict_kwargs_s3: dict = dict()):
        import s3fs

        # save the settings
        self._dict_kwargs_s3 = dict_kwargs_s3

        # open async/sync version of s3fs
        self._as3 = s3fs.S3FileSystem(asynchronous=True, **dict_kwargs_s3)
        self._s3 = s3fs.S3FileSystem(**dict_kwargs_s3)

    def exists(self, path_src: str, **kwargs):
        """# 2023-01-08 23:05:40
        return the list of keys
        """
        return self._s3.exists(path_src, **kwargs)

    def rm(self, path_src: str, flag_recursive: bool = True, **kwargs):
        """# 2023-01-08 23:05:40
        return the list of keys
        """
        return self._s3.rm(path_src, recursive=flag_recursive, **kwargs)  # delete files

    def glob(self, path_src: str, flag_recursive: bool = True, **kwargs):
        """# 2023-01-08 23:05:40
        return the list of keys
        """
        return list(
            "s3://" + e for e in self._s3.glob(path_src, **kwargs)
        )  # 's3://' prefix should be added

    def mkdir(self, path_src: str, **kwargs):
        """# 2023-01-08 23:05:40
        return the list of keys
        """
        # use default 'exist_ok' value
        if "exist_ok" not in kwargs:
            kwargs["exist_ok"] = True
        return self._s3.makedirs(path_src, **kwargs)

    def mv(self, path_src: str, path_dest: str, flag_recursive: bool = True, **kwargs):
        """# 2023-01-08 23:05:40
        return the list of keys
        """
        if not self._s3.exists(
            path_dest, **kwargs
        ):  # avoid overwriting of the existing file
            return self._s3.mv(path_src, path_dest, recursive=flag_recursive, **kwargs)
        else:
            return "destionation file already exists, exiting"

    def cp(self, path_src: str, path_dest: str, flag_recursive: bool = True, **kwargs):
        """# 2023-01-08 23:05:40
        return the list of keys
        """
        if is_s3_url(path_src) and is_s3_url(path_dest):  # copy from s3 to s3
            return self._s3.copy(
                path_src, path_dest, recursive=flag_recursive, **kwargs
            )
        elif is_s3_url(path_src):  # copy from s3 to local
            return self._s3.get(path_src, path_dest, recursive=flag_recursive, **kwargs)
        elif is_s3_url(path_dest):  # copy from local to s3
            return self._s3.put(path_src, path_dest, recursive=flag_recursive, **kwargs)

    def isdir(self, path_src: str, **kwargs):
        """# 2023-01-08 23:05:40
        return the list of keys
        """
        return self._s3.isdir(path_src)

    def get_zarr_metadata(self, path_src: str, **kwargs):
        """# 2023-01-08 23:05:40
        return the list of keys
        ❤️ test
        """
        return dict(zarr.open(path_src).attrs)


class ZarrObject:
    """# 2023-09-24 14:51:46
    A class for hosting Zarr object in a spawned, managed process for accessing remote objects in forked processes
    API functions calls mimic those of a zarr object for seamless replacement of a zarr object.

    path_folder_zarr : str # a path to a (remote) zarr object
    mode : str = 'r' # mode

    path_process_synchronizer : Union[ str, None ] = None # path to the process synchronizer. if None is given, does not use any synchronizer
    """

    def __init__(
        self,
        path_folder_zarr: str,
        mode: str = "r",
        shape: tuple = None,
        chunks: tuple = None,
        dtype=np.int32,
        fill_value=0,
        path_process_synchronizer: Union[str, None] = None,
    ):
        """# 2023-09-24 14:50:36"""
        # set attributes
        self.is_zarr_server = True  # for back-ward compatibility
        self._mode = mode
        self._path_folder_zarr = path_folder_zarr
        self._path_process_synchronizer = path_process_synchronizer

        # open a zarr object
        if mode != "r":  # create a new zarr object, if not available
            if (
                shape is None or chunks is None
            ):  # if one of the arguments for opening zarr array is invalid, open zarr group instead
                za = zarr.open(path_folder_zarr, mode)
            else:  # open zarr array
                za = zarr.open(
                    path_folder_zarr,
                    mode,
                    shape=shape,
                    chunks=chunks,
                    dtype=dtype,
                    fill_value=fill_value,
                )
        else:  # use existing zarr object
            za = zarr.open(path_folder_zarr, mode)
        self._za = za  # set the zarr object as an attribute
        # retrieve attributes of a zarr array
        if hasattr(za, "shape"):  # if zarr object is an array
            self.shape, self.chunks, self.dtype, self.fill_value = (
                self._za.shape,
                self._za.chunks,
                self._za.dtype,
                self._za.fill_value,
            )
        else:  # if zarr object is a group
            self.shape, self.chunks, self.dtype, self.fill_value = (
                None,
                None,
                None,
                None,
            )

    @property
    def flag_hosted(self):
        """# 2023-09-24 14:58:37
        indicates that current object is hosted in a spawned process
        """
        return True

    @property
    def path_folder(self):
        """# 2023-04-19 17:33:21"""
        return self._path_folder_zarr

    def __repr__(self):
        """# 2023-04-20 01:06:16"""
        return f"<HostedZarr of {path_folder_zarr}>"

    @property
    def path_process_synchronizer(self):
        """# 2022-12-07 00:19:29
        return a path of the folder used for process synchronization
        """
        return self._path_process_synchronizer

    def open(
        self,
        path_folder_zarr,
        mode="r",
        shape=None,
        chunks=None,
        dtype=np.int32,
        fill_value=0,
        path_process_synchronizer: Union[str, None] = None,
        reload: bool = False,
    ):
        """# 2023-04-20 02:08:57
        open a new zarr in a ZarrServer object

        reload : bool = False # if True, reload the zarr object even if the 'path_folder' and 'mode' are identical to the currently opened Zarr object. (useful when folder has been updated using the external methods.)
        """
        # if the zarr object is already opened in the same mode, exit, unless 'reload' flag has been set to True.
        if not reload and path_folder_zarr == self.path_folder and self._mode == mode:
            return

        # open a zarr object
        if mode != "r":  # create a new zarr object
            if (
                shape is None or chunks is None
            ):  # if one of the arguments for opening zarr array is invalid, open zarr group instead
                za = zarr.open(path_folder_zarr, mode)
            else:  # open zarr array
                za = zarr.open(
                    path_folder_zarr,
                    mode,
                    shape=shape,
                    chunks=chunks,
                    dtype=dtype,
                    fill_value=fill_value,
                )
        else:  # use existing zarr object
            za = zarr.open(path_folder_zarr, mode)
        self._za = za  # set the zarr object as an attribute
        # retrieve attributes of a zarr array
        if hasattr(za, "shape"):  # if zarr object is an array
            self.shape, self.chunks, self.dtype, self.fill_value = (
                self._za.shape,
                self._za.chunks,
                self._za.dtype,
                self._za.fill_value,
            )
        else:  # if zarr object is a group
            self.shape, self.chunks, self.dtype, self.fill_value = (
                None,
                None,
                None,
                None,
            )
        # update the attributes
        self._path_folder_zarr = path_folder_zarr
        self._mode = mode

    def get_attrs(self, *keys):
        """# 2023-04-19 15:00:04
        get an attribute of the currently opened zarr object using the list of key values
        """
        set_keys = set(self._za.attrs)  # retrieve a list of keys
        return dict(
            (key, self._za.attrs[key]) for key in keys if key in set_keys
        )  # return a subset of metadata using the list of key values given as 'args'

    def get_attr(self, key):
        """# 2023-04-20 01:08:59
        a wrapper of 'get_attrs' for a single key value
        """
        dict_attrs = self.get_attrs(key)  # retrieve the attributes
        if key not in dict_attrs:
            raise KeyError(
                f"attribute {key} does not exist in the zarr object."
            )  # raise a key error if the key does not exist
        return dict_attrs[key]

    def set_attrs(self, **kwargs):
        """# 2023-04-19 15:00:00
        update the attributes of the currently opened zarr object using the keyworded arguments
        """
        # update the metadata
        for key in kwargs:
            self._za.attrs[key] = kwargs[key]

    def get_coordinate_selection(self, *args, **kwargs):
        """# 2022-12-05 22:55:58
        a (possibly) fork-safe wrapper of the 'get_coordinate_selection' zarr operation using a spawned process.
        """
        return self._za.get_coordinate_selection(*args, **kwargs)

    def get_basic_selection(self, *args, **kwargs):
        """# 2022-12-05 22:55:58
        a (possibly) fork-safe wrapper of the 'get_basic_selection' zarr operation using a spawned process.
        """
        return self._za.get_basic_selection(*args, **kwargs)

    def get_orthogonal_selection(self, *args, **kwargs):
        """# 2022-12-05 22:55:58
        a (possibly) fork-safe wrapper of the 'get_orthogonal_selection' zarr operation using a spawned process.
        """
        return self._za.get_orthogonal_selection(*args, **kwargs)

    def get_mask_selection(self, *args, **kwargs):
        """# 2022-12-05 22:55:58
        a (possibly) fork-safe wrapper of the 'get_mask_selection' zarr operation using a spawned process.
        """
        return self._za.get_mask_selection(*args, **kwargs)

    def set_coordinate_selection(self, *args, **kwargs):
        """# 2022-12-05 22:55:58
        a (possibly) fork-safe wrapper of the 'set_coordinate_selection' zarr operation using a spawned process.
        """
        return self._za.set_coordinate_selection(*args, **kwargs)

    def set_basic_selection(self, *args, **kwargs):
        """# 2022-12-05 22:55:58
        a (possibly) fork-safe wrapper of the 'set_basic_selection' zarr operation using a spawned process.
        """
        return self._za.set_basic_selection(*args, **kwargs)

    def set_orthogonal_selection(self, *args, **kwargs):
        """# 2022-12-05 22:55:58
        a (possibly) fork-safe wrapper of the 'set_orthogonal_selection' zarr operation using a spawned process.
        """
        return self._za.set_orthogonal_selection(*args, **kwargs)

    def set_mask_selection(self, *args, **kwargs):
        """# 2022-12-05 22:55:58
        a (possibly) fork-safe wrapper of the 'set_mask_selection' zarr operation using a spawned process.
        """
        return self._za.set_mask_selection(*args, **kwargs)

    def resize(self, *args, **kwargs):
        """# 2022-12-05 22:55:58
        a (possibly) fork-safe wrapper of the 'resize' zarr operation using a spawned process.
        """
        return self._za.resize(*args, **kwargs)

    def __getitem__(self, args):
        """# 2022-12-05 22:55:58
        a (possibly) fork-safe wrapper of the '__getitem__' zarr operation using a spawned process.
        """
        return self._za.__getitem__(args)

    def __setitem__(self, args, values):
        """# 2022-12-05 22:55:58
        a (possibly) fork-safe wrapper of the '__setitem__' zarr operation using a spawned process.
        """
        return self._za.__setitem__(args, values)


# register the manager
ManagerFileSystem.register("FileSystemOperator", FileSystemOperator)
ManagerFileSystem.register("ZarrObject", ZarrObject)
