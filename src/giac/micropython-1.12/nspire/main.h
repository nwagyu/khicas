#ifndef MICROPY_MAIN_H
#define MICROPY_MAIN_H
int do_file(const char *file);
char * micropy_init();
int micropy_eval(const char * line);
void  mp_deinit();
void  mp_stack_ctrl_init();

#endif
