# This file was automatically generated by SWIG (http://www.swig.org).
# Version 1.3.40
#
# Do not make changes to this file unless you know what you are doing--modify
# the SWIG interface file instead.
# This file is compatible with both classic and new-style classes.

from sys import version_info
if version_info >= (2,6,0):
    def swig_import_helper():
        from os.path import dirname
        import imp
        fp = None
        try:
            fp, pathname, description = imp.find_module('_mesh_silhouette_cpp', [dirname(__file__)])
        except ImportError:
            import _mesh_silhouette_cpp
            return _mesh_silhouette_cpp
        if fp is not None:
            try:
                _mod = imp.load_module('_mesh_silhouette_cpp', fp, pathname, description)
            finally:
                fp.close()
            return _mod
    _mesh_silhouette_cpp = swig_import_helper()
    del swig_import_helper
else:
    import _mesh_silhouette_cpp
del version_info
try:
    _swig_property = property
except NameError:
    pass # Python < 2.2 doesn't have 'property'.
def _swig_setattr_nondynamic(self,class_type,name,value,static=1):
    if (name == "thisown"): return self.this.own(value)
    if (name == "this"):
        if type(value).__name__ == 'SwigPyObject':
            self.__dict__[name] = value
            return
    method = class_type.__swig_setmethods__.get(name,None)
    if method: return method(self,value)
    if (not static) or hasattr(self,name):
        self.__dict__[name] = value
    else:
        raise AttributeError("You cannot add attributes to %s" % self)

def _swig_setattr(self,class_type,name,value):
    return _swig_setattr_nondynamic(self,class_type,name,value,0)

def _swig_getattr(self,class_type,name):
    if (name == "thisown"): return self.this.own()
    method = class_type.__swig_getmethods__.get(name,None)
    if method: return method(self)
    raise AttributeError(name)

def _swig_repr(self):
    try: strthis = "proxy of " + self.this.__repr__()
    except: strthis = ""
    return "<%s.%s; %s >" % (self.__class__.__module__, self.__class__.__name__, strthis,)

try:
    _object = object
    _newclass = 1
except AttributeError:
    class _object : pass
    _newclass = 0


class SwigPyIterator(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, SwigPyIterator, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, SwigPyIterator, name)
    def __init__(self, *args, **kwargs): raise AttributeError("No constructor defined - class is abstract")
    __repr__ = _swig_repr
    __swig_destroy__ = _mesh_silhouette_cpp.delete_SwigPyIterator
    __del__ = lambda self : None;
    def value(self): return _mesh_silhouette_cpp.SwigPyIterator_value(self)
    def incr(self, n = 1): return _mesh_silhouette_cpp.SwigPyIterator_incr(self, n)
    def decr(self, n = 1): return _mesh_silhouette_cpp.SwigPyIterator_decr(self, n)
    def distance(self, *args): return _mesh_silhouette_cpp.SwigPyIterator_distance(self, *args)
    def equal(self, *args): return _mesh_silhouette_cpp.SwigPyIterator_equal(self, *args)
    def copy(self): return _mesh_silhouette_cpp.SwigPyIterator_copy(self)
    def next(self): return _mesh_silhouette_cpp.SwigPyIterator_next(self)
    def __next__(self): return _mesh_silhouette_cpp.SwigPyIterator___next__(self)
    def previous(self): return _mesh_silhouette_cpp.SwigPyIterator_previous(self)
    def advance(self, *args): return _mesh_silhouette_cpp.SwigPyIterator_advance(self, *args)
    def __eq__(self, *args): return _mesh_silhouette_cpp.SwigPyIterator___eq__(self, *args)
    def __ne__(self, *args): return _mesh_silhouette_cpp.SwigPyIterator___ne__(self, *args)
    def __iadd__(self, *args): return _mesh_silhouette_cpp.SwigPyIterator___iadd__(self, *args)
    def __isub__(self, *args): return _mesh_silhouette_cpp.SwigPyIterator___isub__(self, *args)
    def __add__(self, *args): return _mesh_silhouette_cpp.SwigPyIterator___add__(self, *args)
    def __sub__(self, *args): return _mesh_silhouette_cpp.SwigPyIterator___sub__(self, *args)
    def __iter__(self): return self
SwigPyIterator_swigregister = _mesh_silhouette_cpp.SwigPyIterator_swigregister
SwigPyIterator_swigregister(SwigPyIterator)

class IntVector(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, IntVector, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, IntVector, name)
    __repr__ = _swig_repr
    def iterator(self): return _mesh_silhouette_cpp.IntVector_iterator(self)
    def __iter__(self): return self.iterator()
    def __nonzero__(self): return _mesh_silhouette_cpp.IntVector___nonzero__(self)
    def __bool__(self): return _mesh_silhouette_cpp.IntVector___bool__(self)
    def __len__(self): return _mesh_silhouette_cpp.IntVector___len__(self)
    def pop(self): return _mesh_silhouette_cpp.IntVector_pop(self)
    def __getslice__(self, *args): return _mesh_silhouette_cpp.IntVector___getslice__(self, *args)
    def __setslice__(self, *args): return _mesh_silhouette_cpp.IntVector___setslice__(self, *args)
    def __delslice__(self, *args): return _mesh_silhouette_cpp.IntVector___delslice__(self, *args)
    def __delitem__(self, *args): return _mesh_silhouette_cpp.IntVector___delitem__(self, *args)
    def __getitem__(self, *args): return _mesh_silhouette_cpp.IntVector___getitem__(self, *args)
    def __setitem__(self, *args): return _mesh_silhouette_cpp.IntVector___setitem__(self, *args)
    def append(self, *args): return _mesh_silhouette_cpp.IntVector_append(self, *args)
    def empty(self): return _mesh_silhouette_cpp.IntVector_empty(self)
    def size(self): return _mesh_silhouette_cpp.IntVector_size(self)
    def clear(self): return _mesh_silhouette_cpp.IntVector_clear(self)
    def swap(self, *args): return _mesh_silhouette_cpp.IntVector_swap(self, *args)
    def get_allocator(self): return _mesh_silhouette_cpp.IntVector_get_allocator(self)
    def begin(self): return _mesh_silhouette_cpp.IntVector_begin(self)
    def end(self): return _mesh_silhouette_cpp.IntVector_end(self)
    def rbegin(self): return _mesh_silhouette_cpp.IntVector_rbegin(self)
    def rend(self): return _mesh_silhouette_cpp.IntVector_rend(self)
    def pop_back(self): return _mesh_silhouette_cpp.IntVector_pop_back(self)
    def erase(self, *args): return _mesh_silhouette_cpp.IntVector_erase(self, *args)
    def __init__(self, *args): 
        this = _mesh_silhouette_cpp.new_IntVector(*args)
        try: self.this.append(this)
        except: self.this = this
    def push_back(self, *args): return _mesh_silhouette_cpp.IntVector_push_back(self, *args)
    def front(self): return _mesh_silhouette_cpp.IntVector_front(self)
    def back(self): return _mesh_silhouette_cpp.IntVector_back(self)
    def assign(self, *args): return _mesh_silhouette_cpp.IntVector_assign(self, *args)
    def resize(self, *args): return _mesh_silhouette_cpp.IntVector_resize(self, *args)
    def insert(self, *args): return _mesh_silhouette_cpp.IntVector_insert(self, *args)
    def reserve(self, *args): return _mesh_silhouette_cpp.IntVector_reserve(self, *args)
    def capacity(self): return _mesh_silhouette_cpp.IntVector_capacity(self)
    __swig_destroy__ = _mesh_silhouette_cpp.delete_IntVector
    __del__ = lambda self : None;
IntVector_swigregister = _mesh_silhouette_cpp.IntVector_swigregister
IntVector_swigregister(IntVector)

class DoubleVector(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, DoubleVector, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, DoubleVector, name)
    __repr__ = _swig_repr
    def iterator(self): return _mesh_silhouette_cpp.DoubleVector_iterator(self)
    def __iter__(self): return self.iterator()
    def __nonzero__(self): return _mesh_silhouette_cpp.DoubleVector___nonzero__(self)
    def __bool__(self): return _mesh_silhouette_cpp.DoubleVector___bool__(self)
    def __len__(self): return _mesh_silhouette_cpp.DoubleVector___len__(self)
    def pop(self): return _mesh_silhouette_cpp.DoubleVector_pop(self)
    def __getslice__(self, *args): return _mesh_silhouette_cpp.DoubleVector___getslice__(self, *args)
    def __setslice__(self, *args): return _mesh_silhouette_cpp.DoubleVector___setslice__(self, *args)
    def __delslice__(self, *args): return _mesh_silhouette_cpp.DoubleVector___delslice__(self, *args)
    def __delitem__(self, *args): return _mesh_silhouette_cpp.DoubleVector___delitem__(self, *args)
    def __getitem__(self, *args): return _mesh_silhouette_cpp.DoubleVector___getitem__(self, *args)
    def __setitem__(self, *args): return _mesh_silhouette_cpp.DoubleVector___setitem__(self, *args)
    def append(self, *args): return _mesh_silhouette_cpp.DoubleVector_append(self, *args)
    def empty(self): return _mesh_silhouette_cpp.DoubleVector_empty(self)
    def size(self): return _mesh_silhouette_cpp.DoubleVector_size(self)
    def clear(self): return _mesh_silhouette_cpp.DoubleVector_clear(self)
    def swap(self, *args): return _mesh_silhouette_cpp.DoubleVector_swap(self, *args)
    def get_allocator(self): return _mesh_silhouette_cpp.DoubleVector_get_allocator(self)
    def begin(self): return _mesh_silhouette_cpp.DoubleVector_begin(self)
    def end(self): return _mesh_silhouette_cpp.DoubleVector_end(self)
    def rbegin(self): return _mesh_silhouette_cpp.DoubleVector_rbegin(self)
    def rend(self): return _mesh_silhouette_cpp.DoubleVector_rend(self)
    def pop_back(self): return _mesh_silhouette_cpp.DoubleVector_pop_back(self)
    def erase(self, *args): return _mesh_silhouette_cpp.DoubleVector_erase(self, *args)
    def __init__(self, *args): 
        this = _mesh_silhouette_cpp.new_DoubleVector(*args)
        try: self.this.append(this)
        except: self.this = this
    def push_back(self, *args): return _mesh_silhouette_cpp.DoubleVector_push_back(self, *args)
    def front(self): return _mesh_silhouette_cpp.DoubleVector_front(self)
    def back(self): return _mesh_silhouette_cpp.DoubleVector_back(self)
    def assign(self, *args): return _mesh_silhouette_cpp.DoubleVector_assign(self, *args)
    def resize(self, *args): return _mesh_silhouette_cpp.DoubleVector_resize(self, *args)
    def insert(self, *args): return _mesh_silhouette_cpp.DoubleVector_insert(self, *args)
    def reserve(self, *args): return _mesh_silhouette_cpp.DoubleVector_reserve(self, *args)
    def capacity(self): return _mesh_silhouette_cpp.DoubleVector_capacity(self)
    __swig_destroy__ = _mesh_silhouette_cpp.delete_DoubleVector
    __del__ = lambda self : None;
DoubleVector_swigregister = _mesh_silhouette_cpp.DoubleVector_swigregister
DoubleVector_swigregister(DoubleVector)

class MeshProject(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, MeshProject, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, MeshProject, name)
    __repr__ = _swig_repr
    def __init__(self, *args): 
        this = _mesh_silhouette_cpp.new_MeshProject(*args)
        try: self.this.append(this)
        except: self.this = this
    __swig_destroy__ = _mesh_silhouette_cpp.delete_MeshProject
    __del__ = lambda self : None;
    def update_vertices(self, *args): return _mesh_silhouette_cpp.MeshProject_update_vertices(self, *args)
    def normals_face(self): return _mesh_silhouette_cpp.MeshProject_normals_face(self)
    def unique_edges(self): return _mesh_silhouette_cpp.MeshProject_unique_edges(self)
    def silhouette_edges(self, *args): return _mesh_silhouette_cpp.MeshProject_silhouette_edges(self, *args)
    def silhouette_edges_lnoe(self, *args): return _mesh_silhouette_cpp.MeshProject_silhouette_edges_lnoe(self, *args)
    def edge_silhouette_num(self): return _mesh_silhouette_cpp.MeshProject_edge_silhouette_num(self)
    def another_test(self): return _mesh_silhouette_cpp.MeshProject_another_test(self)
    def test1(self): return _mesh_silhouette_cpp.MeshProject_test1(self)
    def silhouette_out(self, *args): return _mesh_silhouette_cpp.MeshProject_silhouette_out(self, *args)
MeshProject_swigregister = _mesh_silhouette_cpp.MeshProject_swigregister
MeshProject_swigregister(MeshProject)


