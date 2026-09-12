/*      My image structure has an INFO structure in it and a 2D pixel
        array called DATA. In INFO we have NC, number of columns, and
        NR, the number of rows. Rows are the first index.               */

hough (x, t, error_code)
        struct image *x;
        int t, *error_code;
{
        struct image *z;
        int center_x, center_y, r, omega, i, j, rmax;
        double conv;
        double sin(), cos(), sqrt();

        *error_code = 0;
        conv = 3.1415926535/180.0;
        center_x = x->info->nc/2;       center_y = x->info->nr/2;
        rmax = 
         (int)(sqrt((double)(x->info->nc*x->info->nc+x->info->nr*x->info->nr))/2.0);

/* Create an image for the Hough space - choose your own sampling */
        z = newimage (180, 2*rmax+1, error_code);
        if (*error_code) return;

        for (r = 0; r < 2 * rmax+1; r++)
           for (omega = 0; omega < 180; omega++)
                z->im[r][omega] = 0;

        for (j = 0; j < x->info->nr; j++)
           for (i = 0; i < x->info->nc; i++)
                if (x->data[i][j] <= t)
                   for (omega = 0; omega < 180; ++omega) {
                        r = (i - center_y) * sin((double)(omega*conv)) 
                           + (j - center_x) * cos((double)(omega*conv));
                        if (r>0)r++;
                        z->im[omega][rmax+r] += 1;
                   }
        writeimage (z, "hough.alpha", error_code);
}
