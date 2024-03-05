#include <vector>

typedef std::vector<int> vec1i;
typedef std::vector<double> vec1d;
typedef std::vector<std::vector<int>> vec2i;
typedef std::vector<std::vector<double>> vec2d;
typedef std::vector<std::vector<std::vector<double>>> vec3d;

namespace array_space{

  enum array_data_type = {integer_array,double_array};

  class state_array{
    public:
      state_array(array_data_type, int);
      state_array(array_data_type, int, int, int, int);
  };
}
