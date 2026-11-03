// speed test for printf and std::cout
// $ g++ test/main.cpp -o test/main.out

#include <cstdint>
#include <iostream>
// #include <windows.h>

using namespace std;

long getExecutionTime(void (*func)());
void IfStatement();
void SwitchStatement();
void NoStatement();

int val = 5;
string value = "";

static uint64_t GetTickCountMs()
{
    struct timespec ts;

    clock_gettime(CLOCK_MONOTONIC, &ts);

    return (uint64_t)(ts.tv_nsec / 1000000) + ((uint64_t)ts.tv_sec * 1000ull);
}

long GetExecutionTime(void (*func)()){
    long cTick, nTick;
    cTick = GetTickCountMs();
    (*func)();

    nTick = GetTickCountMs();
    return nTick - cTick;
}

void TestOne(){
    for(int i = 0; i < 100; i++){
        printf("This is the value: %d\n", val);
    }
}

void TestTwo(){
    for(int i = 0; i < 100; i++){
        std::cout << "This is the value: " << val << "\n";
    }
}

void TestThree(){
    //
}

int main()
{
    string testName[] = {"printf", "std::cout"};

    long to[101], tt[101], tth[101];          // [0..99]: runs; [100]: average
    to[100] = 0;
    tt[100] = 0;
    tth[100] = 0;
    std::cout << "Starting test...\n\n";

    for(int i = 0; i < 100; i++)
    {
        to[i] = GetExecutionTime(TestOne);
        tt[i] = GetExecutionTime(TestTwo);
        to[100] += to[i];
        tt[100] += tt[i];
        //tth[100] += tth[i];
        std::cout << "Test is " << i << "% complete...\n\n";
    }

    to[100] /= 100;
    tt[100] /= 100;
    tth[100] /= 100;
    std::cout << "Each function loops 1,000,000 determining the value of an integer and storing it in a string. Each function is executed 100 times to determine the averages.\n\n";
    std::cout << "Function\t\t\tAverage Execution Time (ms)\n";
    std::cout << "-------------------------------------------\n";
    std::cout << "Test One: " << testName[0] << "\t\t" << to[100] << " ms.\n";
    std::cout << "Test Two: " << testName[1] << "\t\t" << tt[100] << " ms.\n\n";
    //std::cout << "Test Three" << testName[2] << "\t\t" << tth[100] << " ms.\n\n";

    system("PAUSE");
    return 0;
}
