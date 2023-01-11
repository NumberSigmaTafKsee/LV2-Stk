# This file was automatically generated by SWIG (http://www.swig.org).
# Version 4.1.0
#
# Do not make changes to this file unless you know what you are doing--modify
# the SWIG interface file instead.

from sys import version_info as _swig_python_version_info
# Import the low-level C/C++ module
if __package__ or "." in __name__:
    from . import _sigpack
else:
    import _sigpack

try:
    import builtins as __builtin__
except ImportError:
    import __builtin__

def _swig_repr(self):
    try:
        strthis = "proxy of " + self.this.__repr__()
    except __builtin__.Exception:
        strthis = ""
    return "<%s.%s; %s >" % (self.__class__.__module__, self.__class__.__name__, strthis,)


def _swig_setattr_nondynamic_instance_variable(set):
    def set_instance_attr(self, name, value):
        if name == "thisown":
            self.this.own(value)
        elif name == "this":
            set(self, name, value)
        elif hasattr(self, name) and isinstance(getattr(type(self), name), property):
            set(self, name, value)
        else:
            raise AttributeError("You cannot add instance attributes to %s" % self)
    return set_instance_attr


def _swig_setattr_nondynamic_class_variable(set):
    def set_class_attr(cls, name, value):
        if hasattr(cls, name) and not isinstance(getattr(cls, name), property):
            set(cls, name, value)
        else:
            raise AttributeError("You cannot add class attributes to %s" % cls)
    return set_class_attr


def _swig_add_metaclass(metaclass):
    """Class decorator for adding a metaclass to a SWIG wrapped class - a slimmed down version of six.add_metaclass"""
    def wrapper(cls):
        return metaclass(cls.__name__, cls.__bases__, cls.__dict__.copy())
    return wrapper


class _SwigNonDynamicMeta(type):
    """Meta class to enforce nondynamic attributes (no new attributes) for a class"""
    __setattr__ = _swig_setattr_nondynamic_class_variable(type.__setattr__)



def sinc(*args):
    return _sigpack.sinc(*args)

def besseli0(x):
    return _sigpack.besseli0(x)

def angle(*args):
    return _sigpack.angle(*args)

def unwrap(x):
    return _sigpack.unwrap(x)

def timevec(N, Fs):
    return _sigpack.timevec(N, Fs)

def sp_version():
    return _sigpack.sp_version()
class FFTW(object):
    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc="The membership flag")
    __repr__ = _swig_repr

    def __init__(self, *args):
        _sigpack.FFTW_swiginit(self, _sigpack.new_FFTW(*args))
    __swig_destroy__ = _sigpack.delete_FFTW

    def fft_cx(self, *args):
        return _sigpack.FFTW_fft_cx(self, *args)

    def ifft_cx(self, *args):
        return _sigpack.FFTW_ifft_cx(self, *args)

    def fft(self, *args):
        return _sigpack.FFTW_fft(self, *args)

    def ifft(self, *args):
        return _sigpack.FFTW_ifft(self, *args)

    def fft2(self, *args):
        return _sigpack.FFTW_fft2(self, *args)

    def ifft2(self, *args):
        return _sigpack.FFTW_ifft2(self, *args)

    def import_wisdom_string(self, wisd):
        return _sigpack.FFTW_import_wisdom_string(self, wisd)

    def import_wisdom_file(self, fname):
        return _sigpack.FFTW_import_wisdom_file(self, fname)

    def export_wisdom_fft(self, fname):
        return _sigpack.FFTW_export_wisdom_fft(self, fname)

    def export_wisdom_ifft(self, fname):
        return _sigpack.FFTW_export_wisdom_ifft(self, fname)

    def export_wisdom_fft_cx(self, fname):
        return _sigpack.FFTW_export_wisdom_fft_cx(self, fname)

    def export_wisdom_ifft_cx(self, fname):
        return _sigpack.FFTW_export_wisdom_ifft_cx(self, fname)

# Register FFTW in _sigpack:
_sigpack.FFTW_swigregister(FFTW)
cvar = _sigpack.cvar
PI = cvar.PI
PI_2 = cvar.PI_2


def cos_win(N, a):
    return _sigpack.cos_win(N, a)

def hamming(N):
    return _sigpack.hamming(N)

def hann(N):
    return _sigpack.hann(N)

def blackman(N):
    return _sigpack.blackman(N)

def blackmanharris(N):
    return _sigpack.blackmanharris(N)

def flattopwin(N):
    return _sigpack.flattopwin(N)

def hanning(N):
    return _sigpack.hanning(N)

def kaiser(N, beta):
    return _sigpack.kaiser(N, beta)

def triang(N):
    return _sigpack.triang(N)

def fir1(M, f0):
    return _sigpack.fir1(M, f0)

def fir1_hp(M, f0):
    return _sigpack.fir1_hp(M, f0)

def fir1_bp(M, f0, f1):
    return _sigpack.fir1_bp(M, f0, f1)

def fir1_bs(M, f0, f1):
    return _sigpack.fir1_bs(M, f0, f1)

def fd_filter(M, fd):
    return _sigpack.fd_filter(M, fd)

def freq(b, a, K=512):
    return _sigpack.freq(b, a, K)

def freqz(b, a, K=512):
    return _sigpack.freqz(b, a, K)

def phasez(b, a, K=512):
    return _sigpack.phasez(b, a, K)
class gplot(object):
    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc="The membership flag")
    __repr__ = _swig_repr

    def __init__(self):
        _sigpack.gplot_swiginit(self, _sigpack.new_gplot())
    __swig_destroy__ = _sigpack.delete_gplot

    def send2gp(self, cmdstr):
        return _sigpack.gplot_send2gp(self, cmdstr)

    def flush_cmd_buf(self):
        return _sigpack.gplot_flush_cmd_buf(self)

    def draw_now(self):
        return _sigpack.gplot_draw_now(self)

    def figure(self, fig):
        return _sigpack.gplot_figure(self, fig)

    def window(self, *args):
        return _sigpack.gplot_window(self, *args)

    def close_window(self):
        return _sigpack.gplot_close_window(self)

    def grid_on(self):
        return _sigpack.gplot_grid_on(self)

    def grid_off(self):
        return _sigpack.gplot_grid_off(self)

    def xlabel(self, label):
        return _sigpack.gplot_xlabel(self, label)

    def ylabel(self, label):
        return _sigpack.gplot_ylabel(self, label)

    def label(self, x, y, label):
        return _sigpack.gplot_label(self, x, y, label)

    def title(self, name):
        return _sigpack.gplot_title(self, name)

    def xlim(self, xmin, xmax):
        return _sigpack.gplot_xlim(self, xmin, xmax)

    def ylim(self, ymin, ymax):
        return _sigpack.gplot_ylim(self, ymin, ymax)

    def plot_add_mat(self, *args):
        return _sigpack.gplot_plot_add_mat(self, *args)

    def plot_show(self):
        return _sigpack.gplot_plot_show(self)

    def plot_clear(self):
        return _sigpack.gplot_plot_clear(self)

    def set_parula_line(self):
        return _sigpack.gplot_set_parula_line(self)

    def set_jet_line(self):
        return _sigpack.gplot_set_jet_line(self)

    def set_set1_line(self):
        return _sigpack.gplot_set_set1_line(self)

    def set_jet_palette(self):
        return _sigpack.gplot_set_jet_palette(self)

    def set_parula_palette(self):
        return _sigpack.gplot_set_parula_palette(self)

    def set_coolwarm_palette(self):
        return _sigpack.gplot_set_coolwarm_palette(self)

    def set_blackbody_palette(self):
        return _sigpack.gplot_set_blackbody_palette(self)

    def set_output(self, name):
        return _sigpack.gplot_set_output(self, name)

    def reset_term(self):
        return _sigpack.gplot_reset_term(self)

    def set_term(self, ttype):
        return _sigpack.gplot_set_term(self, ttype)

# Register gplot in _sigpack:
_sigpack.gplot_swigregister(gplot)

class parser(object):
    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc="The membership flag")
    __repr__ = _swig_repr

    def __init__(self, fname):
        _sigpack.parser_swiginit(self, _sigpack.new_parser(fname))
    __swig_destroy__ = _sigpack.delete_parser

    def getString(self, key, def_val):
        return _sigpack.parser_getString(self, key, def_val)

    def getCxCol(self, key, def_val):
        return _sigpack.parser_getCxCol(self, key, def_val)

    def getCxRow(self, key, def_val):
        return _sigpack.parser_getCxRow(self, key, def_val)

    def getCxMat(self, key, def_val):
        return _sigpack.parser_getCxMat(self, key, def_val)

# Register parser in _sigpack:
_sigpack.parser_swigregister(parser)

class PNM(object):
    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc="The membership flag")
    __repr__ = _swig_repr
    NOTUSED = _sigpack.PNM_NOTUSED
    PBM_A = _sigpack.PNM_PBM_A
    PGM_A = _sigpack.PNM_PGM_A
    PPM_A = _sigpack.PNM_PPM_A
    PBM_B = _sigpack.PNM_PBM_B
    PGM_B = _sigpack.PNM_PGM_B
    PPM_B = _sigpack.PNM_PPM_B
    type = property(_sigpack.PNM_type_get, _sigpack.PNM_type_set)

    def __init__(self):
        _sigpack.PNM_swiginit(self, _sigpack.new_PNM())
    __swig_destroy__ = _sigpack.delete_PNM

    def clear(self):
        return _sigpack.PNM_clear(self)

    def read_header(self):
        return _sigpack.PNM_read_header(self)

    def write_header(self, _type, _rows, _cols, _maxval, comments):
        return _sigpack.PNM_write_header(self, _type, _rows, _cols, _maxval, comments)

    def write(self, *args):
        return _sigpack.PNM_write(self, *args)

    def get_info(self):
        return _sigpack.PNM_get_info(self)

    def get_rows(self):
        return _sigpack.PNM_get_rows(self)

    def get_cols(self):
        return _sigpack.PNM_get_cols(self)

    def get_maxval(self):
        return _sigpack.PNM_get_maxval(self)

    def read(self, *args):
        return _sigpack.PNM_read(self, *args)

# Register PNM in _sigpack:
_sigpack.PNM_swigregister(PNM)


def eval_fcn(*args):
    return _sigpack.eval_fcn(*args)

def lti2discr(F, W, Qc, dT, A, Q):
    return _sigpack.lti2discr(F, W, Qc, dT, A, Q)
class KF(object):
    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc="The membership flag")
    __repr__ = _swig_repr

    def __init__(self, _N, _M, _L):
        _sigpack.KF_swiginit(self, _sigpack.new_KF(_N, _M, _L))
    __swig_destroy__ = _sigpack.delete_KF

    def clear(self):
        return _sigpack.KF_clear(self)

    def set_state_vec(self, _x):
        return _sigpack.KF_set_state_vec(self, _x)

    def set_trans_mat(self, _A):
        return _sigpack.KF_set_trans_mat(self, _A)

    def set_control_mat(self, _B):
        return _sigpack.KF_set_control_mat(self, _B)

    def set_meas_mat(self, _H):
        return _sigpack.KF_set_meas_mat(self, _H)

    def set_err_cov(self, _P):
        return _sigpack.KF_set_err_cov(self, _P)

    def set_proc_noise(self, _Q):
        return _sigpack.KF_set_proc_noise(self, _Q)

    def set_meas_noise(self, _R):
        return _sigpack.KF_set_meas_noise(self, _R)

    def set_kalman_gain(self, _K):
        return _sigpack.KF_set_kalman_gain(self, _K)

    def set_trans_fcn(self, _f):
        return _sigpack.KF_set_trans_fcn(self, _f)

    def set_meas_fcn(self, _h):
        return _sigpack.KF_set_meas_fcn(self, _h)

    def get_state_vec(self):
        return _sigpack.KF_get_state_vec(self)

    def get_err(self):
        return _sigpack.KF_get_err(self)

    def get_kalman_gain(self):
        return _sigpack.KF_get_kalman_gain(self)

    def get_err_cov(self):
        return _sigpack.KF_get_err_cov(self)

    def predict(self, *args):
        return _sigpack.KF_predict(self, *args)

    def update(self, z):
        return _sigpack.KF_update(self, z)

    def rts_smooth(self, Xf, Pf, Xs, Ps):
        return _sigpack.KF_rts_smooth(self, Xf, Pf, Xs, Ps)

# Register KF in _sigpack:
_sigpack.KF_swigregister(KF)

class EKF(KF):
    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc="The membership flag")
    __repr__ = _swig_repr

    def __init__(self, _N, _M, _L):
        _sigpack.EKF_swiginit(self, _sigpack.new_EKF(_N, _M, _L))

    def set_diff_step(self, _dx):
        return _sigpack.EKF_set_diff_step(self, _dx)

    def set_state_jac(self, _f):
        return _sigpack.EKF_set_state_jac(self, _f)

    def set_meas_jac(self, _h):
        return _sigpack.EKF_set_meas_jac(self, _h)

    def jacobian_diff(self, *args):
        return _sigpack.EKF_jacobian_diff(self, *args)

    def jacobian_analytical(self, *args):
        return _sigpack.EKF_jacobian_analytical(self, *args)

    def predict(self, *args):
        return _sigpack.EKF_predict(self, *args)

    def update(self, z):
        return _sigpack.EKF_update(self, z)

    def rts_smooth(self, Xf, Pf, Xs, Ps):
        return _sigpack.EKF_rts_smooth(self, Xf, Pf, Xs, Ps)
    __swig_destroy__ = _sigpack.delete_EKF

# Register EKF in _sigpack:
_sigpack.EKF_swigregister(EKF)

class UKF(KF):
    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc="The membership flag")
    __repr__ = _swig_repr

    def __init__(self, _N, _M, _L):
        _sigpack.UKF_swiginit(self, _sigpack.new_UKF(_N, _M, _L))

    def set_alpha(self, _a):
        return _sigpack.UKF_set_alpha(self, _a)

    def set_beta(self, _b):
        return _sigpack.UKF_set_beta(self, _b)

    def set_kappa(self, _k):
        return _sigpack.UKF_set_kappa(self, _k)

    def set_lambda(self, _l):
        return _sigpack.UKF_set_lambda(self, _l)

    def update_weights(self):
        return _sigpack.UKF_update_weights(self)

    def update_sigma(self, _x, _P):
        return _sigpack.UKF_update_sigma(self, _x, _P)

    def ut(self, _x, _P, _f):
        return _sigpack.UKF_ut(self, _x, _P, _f)

    def predict(self, *args):
        return _sigpack.UKF_predict(self, *args)

    def update(self, z):
        return _sigpack.UKF_update(self, z)

    def rts_smooth(self, Xf, Pf, Xs, Ps):
        return _sigpack.UKF_rts_smooth(self, Xf, Pf, Xs, Ps)
    __swig_destroy__ = _sigpack.delete_UKF

# Register UKF in _sigpack:
_sigpack.UKF_swigregister(UKF)


