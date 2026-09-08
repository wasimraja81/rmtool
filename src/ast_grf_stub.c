/* ast_grf_stub.c -- no-op implementations of the AST "GRF" graphics-
 * primitive interface (Starlink AST's Plot class calls into these; the
 * library itself does not define them, expecting the linking application
 * to supply a real graphics backend or a stub). reproject_cubes never
 * creates an AST Plot object (only FitsChan/FrameSet/Mapping/Resample),
 * so these are never actually invoked at runtime -- they exist purely to
 * satisfy the linker. Signatures taken from /usr/include/grf.h.
 */

int astGAttr(int attr, double value, double *old_value, int prim) {
   (void)attr; (void)value; (void)old_value; (void)prim;
   return 1;
}

int astGScales(float *alpha, float *beta) {
   if (alpha) *alpha = 1.0f;
   if (beta)  *beta  = 1.0f;
   return 1;
}

int astGBBuf(void) { return 1; }
int astGEBuf(void) { return 1; }
int astGFlush(void) { return 1; }

int astGLine(int n, const float *x, const float *y) {
   (void)n; (void)x; (void)y;
   return 1;
}

int astGMark(int n, const float *x, const float *y, int type) {
   (void)n; (void)x; (void)y; (void)type;
   return 1;
}

int astGQch(float *chv, float *chh) {
   if (chv) *chv = 1.0f;
   if (chh) *chh = 1.0f;
   return 1;
}

int astGText(const char *text, float x, float y, const char *just,
             float upx, float upy) {
   (void)text; (void)x; (void)y; (void)just; (void)upx; (void)upy;
   return 1;
}

int astGTxExt(const char *text, float x, float y, const char *just,
              float upx, float upy, float *xb, float *yb) {
   (void)text; (void)x; (void)y; (void)just; (void)upx; (void)upy;
   if (xb) { xb[0]=xb[1]=xb[2]=xb[3]=x; }
   if (yb) { yb[0]=yb[1]=yb[2]=yb[3]=y; }
   return 1;
}

int astGCap(int cap, int value) {
   (void)cap; (void)value;
   return 0;
}
