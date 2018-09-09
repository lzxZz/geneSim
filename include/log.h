#include <string>
#include <initializer_list>
using std::string;
using std::initializer_list;
class Log
{
public:
    static void log(initializer_list<string>);
};