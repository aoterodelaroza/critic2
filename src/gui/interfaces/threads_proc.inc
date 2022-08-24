! -*-F90-*-

!!! tinycthread procedures !!!

!! tinycwrapper
! void *allocate_mtx();
function allocate_mtx() bind(c,name="allocate_mtx")
  import c_ptr
  type(c_ptr) :: allocate_mtx
end function allocate_mtx
!void deallocate_mtx(void *p);
subroutine deallocate_mtx(p) bind(c,name="deallocate_mtx")
  import c_ptr
  type(c_ptr), value :: p
end subroutine deallocate_mtx
!! tinycthread
! int mtx_init(mtx_t *mtx, int type);
function mtx_init(mtx, type) bind(c,name="mtx_init")
  import c_int, c_ptr
  type(c_ptr), value :: mtx
  integer(c_int), value :: type
  integer(c_int) :: mtx_init
end function mtx_init
! void mtx_destroy(mtx_t *mtx);
subroutine mtx_destroy(mtx) bind(c,name="mtx_destroy")
  import c_ptr
  type(c_ptr), value :: mtx
end subroutine mtx_destroy
!x! int mtx_lock(mtx_t *mtx);
!x! int mtx_timedlock(mtx_t *mtx, const struct timespec *ts);
! int mtx_trylock(mtx_t *mtx);
function mtx_trylock(mtx) bind(c,name="mtx_trylock")
  import c_int, c_ptr
  type(c_ptr), value :: mtx
  integer(c_int) :: mtx_trylock
end function mtx_trylock
! int mtx_unlock(mtx_t *mtx);
function mtx_unlock(mtx) bind(c,name="mtx_unlock")
  import c_ptr, c_int
  type(c_ptr), value :: mtx
  integer(c_int) :: mtx_unlock
end function mtx_unlock
!x! int cnd_init(cnd_t *cond);
!x! void cnd_destroy(cnd_t *cond);
!x! int cnd_signal(cnd_t *cond);
!x! int cnd_broadcast(cnd_t *cond);
!x! int cnd_wait(cnd_t *cond, mtx_t *mtx);
!x! int cnd_timedwait(cnd_t *cond, mtx_t *mtx, const struct timespec *ts);
! int thrd_create(thrd_t *thr, thrd_start_t func, void *arg);
function thrd_create(thr, func, arg) bind(c,name="thrd_create")
  import c_ptr, c_funptr, c_int
  type(c_ptr), value :: thr
  type(c_funptr), value :: func
  type(c_ptr), value :: arg
  integer(c_int) :: thrd_create
end function thrd_create
!x! thrd_t thrd_current(void);
!x! int thrd_detach(thrd_t thr);
!x! int thrd_equal(thrd_t thr0, thrd_t thr1);
!x! TTHREAD_NORETURN void thrd_exit(int res);
! int thrd_join(thrd_t thr, int *res); --- wrapped ---
! int wrap_thrd_join(thrd_t *thr, int *res)
function wrap_thrd_join(thr, res) bind(c,name="wrap_thrd_join")
  import c_int, c_ptr
  type(c_ptr), value :: thr
  integer(c_int) :: res
  integer(c_int) :: wrap_thrd_join
end function wrap_thrd_join
!x! int thrd_sleep(const struct timespec *duration, struct timespec *remaining);
!x! void thrd_yield(void);
!x! int tss_create(tss_t *key, tss_dtor_t dtor);
!x! void tss_delete(tss_t key);
!x! void *tss_get(tss_t key);
!x! int tss_set(tss_t key, void *val);
