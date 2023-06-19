#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <memory.h>
#include "malloc.h"
#include "libs/algorithms/algorithms.h"

#define ARRAY_SIZE(a) sizeof(a) / sizeof(a[0])

#define TIME_TEST(testCode, time) { \
    clock_t start_time = clock(); \
    testCode \
        clock_t end_time = clock(); \
    clock_t sort_time = end_time - start_time; \
    time = (double) sort_time / CLOCKS_PER_SEC; \
}



// функция сортировки
typedef struct SortFunc {
    void (*sort )(int *a, size_t n); // указатель на функцию
    // сортировки
    char name[64];                   // имя сортировки,
    // используемое при выводе
} SortFunc;

// функция генерации
typedef struct GeneratingFunc {
    void (*generate )(int *a, size_t n); // указатель на функцию
    // генерации последоват.
    char name[64];                       // имя генератора,
    // используемое при выводе
} GeneratingFunc;

// функция сортировки
typedef struct nCompSort {
    long long (*nComp )(int *a, size_t n); // указатель на функцию
    // сортировки
    char name[64];                   // имя сортировки,
    // используемое при выводе
} nCompSort;

void checkTime(void (*sortFunc)(int *, size_t),
               void (*generateFunc)(int *, size_t),
               size_t size, char *experimentName) {
    static size_t runCounter = 1;

    // генерация последовательности
    static int innerBuffer[100000];
    generateFunc(innerBuffer, size);
    printf("Run #%zu| ", runCounter++);
    printf("Name: %s\n", experimentName);

    // замер времени
    double time;

    TIME_TEST({ sortFunc(innerBuffer, size); }, time);

    // результаты замера
    printf("Status: ");
    if (isOrdered(innerBuffer, size)) {
        printf("OK! Time : %.3f s.\n", time);

        // запись в файл
        char filename[256];
        sprintf(filename, "./data/%s.csv", experimentName);
        FILE *f = fopen(filename, "a");
        if (f == NULL) {
            printf("FileOpenError %s", filename);
            exit(1);
        }
        fprintf(f, "%zu; %.3f\n", size, time);
        fclose(f);
    } else {
        printf("Wrong!\n");

        // вывод массива, который не смог быть отсортирован
        outputArray_(innerBuffer, size);

        exit(1);
    }
}

void checkNComp(long long (*nComp )(int *a, size_t n),
                void (*generateFunc)(int *, size_t),
                size_t size, char *experimentName, char *name) {
    static size_t runCounter = 1;

    // генерация последовательности
    static int innerBuffer[100000];
    generateFunc(innerBuffer, size);
    printf("Run #%zu| ", runCounter++);
    printf("Name: %s\n", experimentName);

    // замер времени
    long long nComps = nComp(innerBuffer, size);

    // результаты замера
    printf("Status: ");
    if (isOrdered(innerBuffer, size)) {
        printf("OK! Comps: %lld\n", nComps);

        // запись в файл
        char filename[256];
        sprintf(filename, "./data/%s.csv",experimentName);
        FILE *f = fopen(filename, "a");
        if (f == NULL) {
            printf("FileOpenError %s", filename);
            exit(1);
        }
        fprintf(f, "%zu; %lld\n", size, nComps);
        fclose(f);
    } else {
        printf("Wrong!\n");

        // вывод массива, который не смог быть отсортирован
        outputArray_(innerBuffer, size);

        exit(1);
    }
}


void selectionSort(int *a, size_t size) {
    for (int i = 0; i < size - 1; i++) {
        int minPos = i;
        for (int j = i + 1; j < size; j++)
            if (a[j] < a[minPos])
                minPos = j;
        swap(&a[i], &a[minPos], sizeof(int));
    }
}

long long selectionSortN(int *a, size_t size) {
    long long countComp = 0;
    for (int i = 0; i < size - 1 && ++countComp; i++) {
        int minPos = i;
        for (int j = i + 1; j < size && ++countComp; j++)
            if (a[j] < a[minPos] && ++countComp)
                minPos = j;
        swap(&a[i], &a[minPos], sizeof(int));
    }
    return countComp;
}


void insertionSort(int *a, size_t size) {
    for (size_t i = 1; i < size; i++) {
        int t = a[i];
        size_t j = i;
        while (j > 0 && a[j - 1] > t) {
            a[j] = a[j - 1];
            j--;
        }
        a[j] = t;
    }
}

long long insertionSortN(int *a, size_t size) {
    long long countComp = 0;
    for (size_t i = 1; i < size && ++countComp; i++) {
        int t = a[i];
        size_t j = i;
        while (j > 0 && ++countComp && a[j - 1] > t && ++countComp) {
            a[j] = a[j - 1];
            j--;
        }
        a[j] = t;
    }
    return countComp;
}

void bubbleSort(int *a, size_t size) {
    for (size_t i = 0; i < size - 1; i++)
        for (size_t j = size - 1; j > i; j--)
            if (a[j - 1] > a[j])
                swap(&a[j - 1], &a[j], sizeof(int));
}

long long bubbleSortN(int *a, size_t size) {
    long long countComp = 0;
    for (size_t i = 0; i < size - 1 && ++countComp; i++)
        for (size_t j = size - 1; j > i && ++countComp; j--)
            if (a[j - 1] > a[j] && ++countComp)
                swap(&a[j - 1], &a[j], sizeof(int));

    return countComp;
}

void combsort(int *a, const size_t size) {
    size_t step = size;
    int swapped = 1;
    while (step > 1 || swapped) {
        if (step > 1)
            step /= 1.24733;
        swapped = 0;
        for (size_t i = 0, j = i + step; j < size; ++i, ++j)
            if (a[i] > a[j]) {
                swap(&a[i], &a[j], sizeof(int));
                swapped = 1;
            }
    }
}

long long combsortN(int *a, const size_t size) {
    size_t step = size;
    int swapped = 1;
    long long countComp = 0;
    while (step > 1 && ++countComp || swapped && ++countComp) {
        if (step > 1 && ++countComp)
            step /= 1.24733;
        swapped = 0;
        for (size_t i = 0, j = i + step; j < size && ++countComp; ++i, ++j)
            if (a[i] > a[j] && ++countComp) {
                swap(&a[i], &a[j], sizeof(int));
                swapped = 1;
            }
    }
    return countComp;
}

void shellSort(int *a, size_t size) {
    for (size_t step = size / 2; step > 0; step /= 2)
        for (size_t i = step; i < size; i++) {
            int tmp = a[i];
            size_t j;
            for (j = i; j >= step; j -= step) {
                if (tmp < a[j - step])
                    a[j] = a[j - step];
                else
                    break;
            }
            a[j] = tmp;
        }
}

long long shellSortN(int *a, size_t size) {
    long long countComp = 0;
    for (size_t step = size / 2; step > 0 && ++countComp; step /= 2)
        for (size_t i = step; i < size && ++countComp; i++) {
            int tmp = a[i];
            size_t j;
            for (j = i; j >= step && ++countComp; j -= step) {
                if (tmp < a[j - step]&&++countComp)
                    a[j] = a[j - step];
                else
                    break;
            }
            a[j] = tmp;
        }
    return countComp;
}


int cmp_int(const void *a, const void *b) {
    return *(const int *) a - *(const int *) b;
}

int cmp_intReverse(const void *a, const void *b) {
    return *(const int *) b - *(const int *) a;
}

int digit(int n, int k, int N, int M) {
    return (n >> (N * k) & (M - 1));
}

void _radixSort(int *l, int *r, int N) {
    int k = (32 + N - 1) / N;
    int M = 1 << N;
    int sz = r - l;
    int *b = (int *) malloc(sizeof(int) * sz);
    int *c = (int *) malloc(sizeof(int) * M);
    for (int i = 0; i < k; i++) {
        for (int j = 0; j < M; j++)
            c[j] = 0;

        for (int *j = l; j < r; j++)
            c[digit(*j, i, N, M)]++;

        for (int j = 1; j < M; j++)
            c[j] += c[j - 1];

        for (int *j = r - 1; j >= l; j--)
            b[--c[digit(*j, i, N, M)]] = *j;

        int cur = 0;
        for (int *j = l; j < r; j++)
            *j = b[cur++];
    }
    free(b);
    free(c);
}

void radix_Sort(int *a, size_t n) {
    _radixSort(a, a + n, 8);
}
void rad_sort_u(unsigned *begin, unsigned *end, unsigned bit) {
    if (!bit || end < begin + 1)
        return;

    unsigned *ll = begin;
    unsigned *rr = end - 1;
    while (1) {
        while (ll < rr && !(*ll & bit))
            ll++;
        while (ll < rr && (*rr & bit))
            rr--;
        if (ll >= rr)
            break;
        swap(ll, rr, sizeof(int));
    }

    if (!(bit & *ll) && ll < end)
        ll++;

    bit >>= 1;

    rad_sort_u(begin, ll, bit);
    rad_sort_u(ll, end, bit);
}

long long rad_sort_uN(unsigned *begin, unsigned *end, unsigned bit) {
    long long countComp = 0;
    if (!bit && ++countComp || end < begin + 1 && ++countComp)
        return countComp;

    unsigned *ll = begin;
    unsigned *rr = end - 1;
    while (++countComp) {
        while (ll < rr && ++countComp && !(*ll & bit) && ++countComp)
            ll++;
        while (ll < rr && ++countComp && (*rr & bit) && ++countComp)
            rr--;
        if (ll >= rr && ++countComp)
            break;
        swap(ll, rr, sizeof(int));
    }

    if (!(bit & *ll) && ++countComp && ll < end && ++countComp)
        ll++;

    bit >>= 1;

    rad_sort_uN(begin, ll, bit);
    rad_sort_uN(ll, end, bit);

    return countComp;
}

long long radix_sortN(int *a, const size_t size) {
    unsigned *x = (unsigned *) a;

    long long countComp = 0;
    for (size_t i = 0; i < size && ++countComp; i++)
        x[i] ^= INT_MIN;

    countComp += rad_sort_uN(x, x + size, INT_MIN);

    for (size_t i = 0; i < size && ++countComp; i++)
        x[i] ^= INT_MIN;

    return countComp;
}

void merge(const int *a, const size_t n,
           const int *b, const size_t m, int *c) {
    int i = 0, j = 0;
    while (i < n || j < m) {
        if (j == m || i < n && a[i] < b[j]) {
            c[i + j] = a[i];
            i++;
        } else {
            c[i + j] = b[j];
            j++;
        }
    }
}

long long mergeN(const int *a, const size_t n,
                 const int *b, const size_t m, int *c) {
    long long countComp = 0;
    int i = 0, j = 0;
    while (i < n && ++countComp || j < m && ++countComp) {
        if (j == m && ++countComp || i < n && ++countComp && a[i] < b[j] && ++countComp) {
            c[i + j] = a[i];
            i++;
        } else {
            c[i + j] = b[j];
            j++;
        }
    }
    return countComp;
}


void mergeSort_(int *source, size_t l, size_t r, int *buffer) {
    size_t n = r - l;
    if (n <= 1)
        return;

    size_t m = (l + r) / 2;
    mergeSort_(source, l, m, buffer);
    mergeSort_(source, m, r, buffer);

    merge(source + l, m - l, source + m, r - m, buffer);
    memcpy(source + l, buffer, sizeof(int) * n);
}

long long mergeSortN_(int *source, size_t l, size_t r, int *buffer) {
    size_t n = r - l;
    long long countComp = 0;
    if (n <= 1 && ++countComp)
        return countComp;

    size_t m = (l + r) / 2;
    mergeSort_(source, l, m, buffer);
    mergeSort_(source, m, r, buffer);

    countComp = mergeN(source + l, m - l, source + m, r - m, buffer);
    memcpy(source + l, buffer, sizeof(int) * n);

    return countComp;
}


void mergeSort(int *a, size_t n) {
    int *buffer = (int *) malloc(sizeof(int) * n);
    mergeSort_(a, 0, n, buffer);
    free(buffer);
}

long long mergeSortN(int *a, size_t n) {
    int *buffer = (int *) malloc(sizeof(int) * n);
    long long countComp = mergeSortN_(a, 0, n, buffer);
    free(buffer);

    return countComp;
}


void quicksort(int *A, size_t len) {
    if (len < 2)
        return;

    int pivot = A[len / 2];

    int i, j;
    for (i = 0, j = len - 1;; i++, j--) {
        while (A[i] < pivot)
            i++;
        while (A[j] > pivot)
            j--;

        if (i >= j)
            break;

        swap(A + i, A + j, sizeof(int));
    }

    quicksort(A, i);
    quicksort(A + i, len - i);
}

long long quicksortN(int *A, size_t len) {
    long long countComp = 0;
    if (len < 2 && ++countComp)
        return countComp;

    int pivot = A[len / 2];

    int i, j;
    for (i = 0, j = len - 1; ++countComp; i++, j--) {
        while (A[i] < pivot && ++countComp)
            i++;
        while (A[j] > pivot && ++countComp)
            j--;

        if (i >= j && ++countComp)
            break;

        swap(A + i, A + j, sizeof(int));
    }

    countComp += quicksortN(A, i);
    countComp += quicksortN(A + i, len - i);

    return countComp;
}

void generateRandomArray(int *a, size_t size) {
    srand(time(0));
    for (int i = 0; i < size; i++)
        a[i] = 20000 - rand() % 10000;
}

void generateOrderedArray(int *a, size_t size) {
    generateRandomArray(a, size);
    qsort(a, size, sizeof(int), cmp_int);
}

void generateOrderedBackwards(int *a, size_t size) {
    generateRandomArray(a, size);
    qsort(a, size, sizeof(int), cmp_intReverse);
}

void qSort(int *a, size_t size) {
    qsort(a, size, sizeof(int), cmp_int);
}

void heapify(int *arr, int n, int i) {
    int largest = i;
    int l = 2 * i + 1;
    int r = 2 * i + 2;

    if (l < n && arr[l] > arr[largest])
        largest = l;
    if (r < n && arr[r] > arr[largest])
        largest = r;
    if (largest != i) {
        swap(&arr[i], &arr[largest], sizeof(int));
        heapify(arr, n, largest);
    }
}

void heapSort(int *arr, size_t n) {
    for (int i = n / 2 - 1; i >= 0; i--)
        heapify(arr, n, i);
    for (int i = n - 1; i >= 0; i--) {
        swap(&arr[0], &arr[i], sizeof(int));
        heapify(arr, i, 0);
    }
}

long long heapifyN(int *arr, int n, int i) {
    long long countComp = 0;
    int largest = i;
    int l = 2 * i + 1;
    int r = 2 * i + 2;

    if (l < n && ++countComp && arr[l] > arr[largest] && ++countComp)
        largest = l;

    if (r < n && ++countComp && arr[r] > arr[largest] && ++countComp)
        largest = r;

    if (largest != i && ++countComp) {
        swap(&arr[i], &arr[largest], sizeof(int));

        countComp += heapifyN(arr, n, largest);
    }

    return countComp;
}

long long heapSortN(int *arr, size_t n) {
    long long countComp = 0;
    for (int i = n / 2 - 1; i >= 0 && ++countComp; i--)
        countComp += heapifyN(arr, n, i);

    for (int i = n - 1; i >= 0 && ++countComp; i--) {
        swap(&arr[0], &arr[i], sizeof(int));

        countComp += heapifyN(arr, i, 0);
    }
    return countComp;
}

////////////////////////////////////////////////////

void timeExperiment() {
    nCompSort nComps[] = {
            {selectionSortN, "selectionSortN"},
            {insertionSortN, "insertionSortN"},
            {bubbleSortN,    "bubbleSortN"},
            {combsortN,      "combSortN"},
            {shellSortN,     "shellSortN"},
            {radix_sortN,    "radixSortN"},
            {mergeSortN,     "mergeSortN"},
            {heapSortN,      "heapSortN"},
            {quicksortN,     "quickSortN"},
    };

    const unsigned COMPS_N = ARRAY_SIZE(nComps);

    // описание функций генерации
    GeneratingFunc generatingFuncs[] = {
            {generateRandomArray,      "random"},
            {generateOrderedArray,     "ordered"},
            {generateOrderedBackwards, "orderedBackwards"}
    };
    const unsigned CASES_N = ARRAY_SIZE(generatingFuncs);

    // запись статистики в файл
//    for (size_t size = 10000; size <= 100000; size += 10000) {
//        printf(" - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n");
//        printf("Size : %d\n", size);
//        for (int i = 0; i < FUNCS_N; i++) {
//            for (int j = 0; j < CASES_N; j++) {
//                // генерация имени файла
//                static char filename[128];
//                sprintf(filename, "%s_%s_time",
//                        sorts[i].name, generatingFuncs[j].name);
//                checkTime(sorts[i].sort,
//                          generatingFuncs[j].generate,
//                          size, filename, sorts[i].name);
//            }
//        }
//        printf("\n");
//    }

    // запись статистики в файл
    for (size_t size = 10000; size <= 100000; size += 10000) {
        printf(" - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n");
        printf("Size : %d\n", size);
        for (int i = 0; i < COMPS_N; i++) {
            for (int j = 0; j < CASES_N; j++) {
                // генерация имени файла
                static char filename[128];
                sprintf(filename, "%s_%s_time",
                        nComps[i].name, generatingFuncs[j].name);
                checkNComp(nComps[i].nComp,
                           generatingFuncs[j].generate,
                           size, filename, nComps[i].name);
            }
        }
        printf("\n");
    }
}

int main() {
    timeExperiment();

    return 0;
}