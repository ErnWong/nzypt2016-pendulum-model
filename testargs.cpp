#include <cstdio>

int main(int argc, char ** argv)
{
    printf("argc: %d\n", argc);
    for (int i = 0; i < argc; i++)
    {
        printf("argv[i]: %s\n", argv[i]);
    }
}
