! -*-F90-*-

!xx! tinycthread procedures

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
! void *allocate_thrd();
function allocate_thrd() bind(c,name="allocate_thrd")
  import c_ptr
  type(c_ptr) :: allocate_thrd
end function allocate_thrd
!void deallocate_thrd(void *p);
subroutine deallocate_thrd(p) bind(c,name="deallocate_thrd")
  import c_ptr
  type(c_ptr), value :: p
end subroutine deallocate_thrd
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
!xx! int mtx_lock(mtx_t *mtx);
!xx! int mtx_timedlock(mtx_t *mtx, const struct timespec *ts);
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
!xx! int cnd_init(cnd_t *cond);
!xx! void cnd_destroy(cnd_t *cond);
!xx! int cnd_signal(cnd_t *cond);
!xx! int cnd_broadcast(cnd_t *cond);
!xx! int cnd_wait(cnd_t *cond, mtx_t *mtx);
!xx! int cnd_timedwait(cnd_t *cond, mtx_t *mtx, const struct timespec *ts);
! int thrd_create(thrd_t *thr, thrd_start_t func, void *arg);
function thrd_create(thr, func, arg) bind(c,name="thrd_create")
  import c_ptr, c_funptr, c_int
  type(c_ptr), value :: thr
  type(c_funptr), value :: func
  type(c_ptr), value :: arg
  integer(c_int) :: thrd_create
end function thrd_create
!xx! thrd_t thrd_current(void);
!xx! int thrd_detach(thrd_t thr);
!xx! int thrd_equal(thrd_t thr0, thrd_t thr1);
!xx! TTHREAD_NORETURN void thrd_exit(int res);
!xx! int thrd_join(thrd_t thr, int *res);
!xx! int thrd_sleep(const struct timespec *duration, struct timespec *remaining);
!xx! void thrd_yield(void);
!xx! int tss_create(tss_t *key, tss_dtor_t dtor);
!xx! void tss_delete(tss_t key);
!xx! void *tss_get(tss_t key);
!xx! int tss_set(tss_t key, void *val);
