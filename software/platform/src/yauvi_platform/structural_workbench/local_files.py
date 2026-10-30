"""Never trigger implicit hydration while fingerprinting optional local source."""
import sys

class LocalSourceUnavailable(OSError):
    pass

def require_local(path):
    # Darwin sys/stat.h defines SF_DATALESS as 0x40000000. An ordinary read of
    # such a File Provider placeholder can wait indefinitely for hydration.
    info=path.stat()
    if sys.platform=='darwin' and getattr(info,'st_flags',0)&0x40000000:
        raise LocalSourceUnavailable('Local source is offloaded: '+path.name)

def local_bytes(path):
    require_local(path)
    return path.read_bytes()
