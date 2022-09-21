#include "datafiltering.h"
#include "ASR.h"
#include "2d_array.h"
#include "online_ASR.h"
#include "offline_ASR.h"

int main()
{
    //online_ASR();
    offline_ASR(500, 8, 40000);
    return 0;
}


