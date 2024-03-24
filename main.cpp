#include <iostream>

class Oscillazione {
    public:
        Oscillazione(){
            std::cout << "sesso pazzo oscillatorio" << std::endl;
        }
};

int main(int argv, char** argc){
    std::cout << "sesso pazzo" << std::endl;
    Oscillazione* n1 = new Oscillazione;
    return 0;
}