#include <stdio.h>
#include <stdlib.h>
#include "tcp_client.h"

using namespace std;

int main(int argc,char* argv[])
{
        tcp_client tc(argv[1],argv[2], argv[3], argv[4]);
        return 0;
}

