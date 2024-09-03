#include <iostream>

int main(int argc, char* argv[]) {
	
	int N = std::stoi(argv[1]);

	int a = 0;

	for (int i = 0; i < N; ++i) {
        a += 1;
    }

    std::cout << a << std::endl;

}