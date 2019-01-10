#include <fstream>
#include <string>
#include <cstdio>
int main(int arcc, const char* argv[]) {
  std::ifstream in;
  in.open(argv[1]);
  if (!in.is_open()) {
    fprintf(stderr, "Fail to open file: %s\n", argv[1]);
  }
  std::string line;
  double value;
  int i = 0;
  while (std::getline(in, line) && line != "") {
    sscanf(line.c_str(), "%lf", &value);
    if ((i + 1) % 12 == 0) printf("%lf\n", value);
    else printf("%lf ", value);
    ++i;
  }
  return 0;
}
