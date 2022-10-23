import ctypes
mes = ctypes.CDLL('./fib.so');
mes.main.argtypes = (ctypes.c_int, ctypes.POINTER(ctypes.c_char_p))
def measure(a):
    global mes
    length=len(a)
    arr_type = ctypes.c_char_p*length
    valueArr = ctypes.c_double
    result = mes.main(ctypes.c_int(length), arr_type(*a))
    return result
