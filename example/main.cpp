#include <string>
#include <vector>
#include <windows.h>

extern "C" int utf8_main(int, char **);

std::string to_utf8(const std::wstring &ws)
{
    int size = WideCharToMultiByte(65001, 0, ws.c_str(), -1, 0, 0, 0, 0);
    std::vector<char> buffer(size);
    size = WideCharToMultiByte(65001, 0, ws.c_str(), -1, &buffer[0], size, 0, 0);
    return std::string(buffer.begin(), buffer.end());
}

int wmain1(int argc, wchar_t **argv)
{
    std::vector<std::string> args(argc);
    for (int i = 0; i < argc; ++i)
	args[i] = to_utf8(argv[i]);
    std::vector<char*> cargs(argc + 1);
    for (int i = 0; i < argc; ++i)
	cargs[i] = const_cast<char*>(args[i].c_str());
    return utf8_main(argc, &cargs[0]);
}

typedef struct { int newmode; } _startupinfo;
extern "C"
int __wgetmainargs(int *, wchar_t ***, wchar_t ***, int, _startupinfo *);

int main()
{
    int argc;
    wchar_t **argv, **envp;
    _startupinfo si = { 0 };
    __wgetmainargs(&argc, &argv, &envp, 1, &si);
    return wmain1(argc, argv);
}
