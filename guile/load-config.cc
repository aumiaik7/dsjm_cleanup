#include <cstdio>
#include <guile/gh.h>

#define CONFIGFILENAME ".gcolor"
#define DIRECTORYSEPARATOR "/"

#include "Configuration.hh"
#include "driver.hh"

Configuration *configuration;

void read_config_file(void)
{
  configuration = new Configuration();
  char *filename;
  filename = (char *) malloc(strlen(CONFIGFILENAME) + 1);
  if(filename == NULL) return;
  sprintf(filename, "%s",CONFIGFILENAME);

  gh_eval_file(filename);
  free(filename);
}

SCM hello_world()
{
  printf("hello world\n");
  return SCM_EOL;
}

void register_procs(void)
{
  gh_new_procedure("hello-world", hello_world,0,0,0);
}

void inner_main(int argc, char **argv)
{
  register_procs();
  read_config_file();
  printf("... Configuration Loaded\n");
  run(configuration);
}

int main(int argc, char **argv)
{
  printf("gcolor-load-config\n");
  gh_enter(argc,argv,inner_main);

  return 0;
}
