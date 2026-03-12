/*
 * This is used to fail the build tests of someone tries to run them without
 * enabling them in the build first.
 */
#include <stdio.h>

int main() {
    printf("No tests run. Need to enable building tests with meson configure -Dtests=enabled\n");
    return 1;
}
