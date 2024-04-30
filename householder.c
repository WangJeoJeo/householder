#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

#define MAX 1024
#define ROW 4
#define COL 4

// 返回x的绝对值，y的符号
double sign(double x, double y)
{
    return fabs(x) * (y >= 0 ? 1 : -1);
}
// 获取每列元素，u,s,beta
void get_ele(double **src, int n, int row, double u[], double *s, double *beta)
{
    int i, j;
    double s1;
    *beta = *s = s1 = 0.0;
    for (i = n; i < row; i++)
    {
        s1 += src[i][n] * src[i][n];
    }
    s1 = -sqrt(s1) * sign(1, src[n][n]);
    if (s1 == 0.0)
    {
        printf("*****WARRING(get_ele):s equal 0.0, column: %d\n", n + 1);
        return;
    }
    u[0] = src[n][n] - s1;
    if (u[0] == 0.0)
    {
        printf("*****WARRING(get_ele):u[0] equal 0.0, column: %d\n", n + 1);
        return;
    }
    for (i = 1; i < row - n; i++)
    {
        u[i] = src[n + i][n];
    }
    *beta = 1.0 / s1 / u[0];
    *s = s1;
}
// householder变换
void householder(double **src, int row, int col, bool save)
{
    int i, j, k, flag = 1, tmp;
    double s, beta, gama;
    double u[MAX] = {0};

    if (src == NULL)
    {
        printf("*****ERROR(householder):src Matrix is NULL, exit!\n");
        exit(-1);
    }
    for (i = 0; i < col; i++)
    {
        // check elemental first
        flag = 1;
        tmp = 0;
        for (j = i + 1; j < row; j++)
        {
            if (src[j][i] == 0)
            {
                tmp++;
            }
        }
        if (tmp == row)
        {
            flag = 0;
            if (save)
            {
                src[j + 1][i] = 0;
                src[j = 2][i] = 0;
            }
        }
        if (flag)
        {
            get_ele(src, i, col, u, &s, &beta);
            if (s == 0.0)
                continue;
            src[i][i] = s;
            for (j = i + 1; j < col; j++)
            {
                gama = 0.0;
                for (k = i; k < row; k++)
                {
                    gama += u[k - i] * src[k][j];
                }
                if (gama == 0.0)
                    continue;
                gama *= beta;
                for (k = i; k < row; k++)
                {
                    src[k][j] += gama * u[k - i];
                }
            }
        }
        if (save)
        {
            for (j = i + 1; j < row + 1; j++)
            {
                src[j][i] = u[j - i - 1];
            }
            src[j][i] = beta;
        }
        else
        {
            for (j = i + 1; j < row; j++)
            {
                src[j][i] = 0.0;
            }
        }
    }
}

void test1(bool save)
{
    int i, j;
    double **arr = (double **)malloc(sizeof(double *) * (ROW + 2));
    for (i = 0; i < ROW + 2; i++)
    {
        arr[i] = (double *)malloc(sizeof(double) * COL);
        memset(arr[i], 0, sizeof(double) * COL);
    }
    arr[0][0] = 1;
    arr[1][0] = 1;
    arr[2][0] = 1;
    arr[3][0] = 1;

    arr[0][1] = 2;
    arr[1][1] = 0;
    arr[2][1] = 0;
    arr[3][1] = 2;

    arr[0][2] = 0;
    arr[1][2] = 3;
    arr[2][2] = 3;
    arr[3][2] = 0;

    arr[0][3] = 1;
    arr[1][3] = 1;
    arr[2][3] = 2;
    arr[3][3] = 2;

    householder(arr, ROW, COL, save);

    for (i = 0; i < ROW + 2; i++)
    {
        for (j = 0; j < COL; j++)
        {
            printf("%5.2lf  ", arr[i][j]);
        }
        putchar('\n');
    }
}

int main(int argc, char *argv[])
{
    test1(false);
    test1(true);
    return 0;
}