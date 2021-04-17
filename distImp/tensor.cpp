#include <iostream>
#include <vector>

using namespace std;

template<typename T>
class Tensor
{
    public:
        //copy constructor
        Tensor(const Tensor&) = default;
        //copy assignment operator
        Tensor& operator=(Tensor&) = default;
        //move constructor
        //Tensor(const Tensor&&) = default;
        //move assignment operator
        Tensor& operator=(Tensor&&) = default;
        template<typename... Dimensions>
        Tensor(Dimensions... dims)
        {
            unsigned int nelements = helper_constructor(dims...);
            data = new T [nelements];
        }
        
        unsigned int helper_constructor(unsigned int dim)
        {
            return dim;
        }
        template<typename... Dimensions>
        unsigned int helper_constructor(unsigned int dim, Dimensions... dims)
        {
            return dim * helper_constructor(dims...);
        }
        Tensor(std::initializer_list<unsigned int> list)
        {

            unsigned int nelements = 1;
            for( auto elem : list )
            {
                dimensions.push_back(elem);
                nelements = nelements * elem;
            }
            data = new T [nelements];
        }
        ~Tensor()
        {
            if(data != nullptr) delete [] data;
        }
    private:
        T* data = nullptr;
        vector<unsigned int> dimensions;
};


int main()
{
    Tensor<int> A({2,3,4});
    Tensor<int> B(3,4);
    return 0;
}

