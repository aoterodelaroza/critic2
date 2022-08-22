#if defined(_WIN32) || defined(__WIN32__) || defined(__WINDOWS__)
  /* Windows */
  #include <direct.h>
  #define GETCWD _getcwd
#else
  /* Unix */
  #include <unistd.h>
  #define GETCWD getcwd
#endif

int getCurrentWorkDir(char *str, size_t siz){
  if (GETCWD(str, siz) == str)
    return 0;
  else
    return 1;
}
