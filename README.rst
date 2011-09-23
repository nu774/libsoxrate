libsoxrate
==========

What is this?
-------------
Win32 DLL based on rate converter of SoX(http://sox.sourceforge.net/).
Difference from libsox are:
- Only contains rate converter.
- Not a static link lib.
- Different API.
- Using Win32 native threads for parallelizing.
- Eliminated needs of synchronization for global cache in the original code.

However, basically it's nothing but sox rate converter.
