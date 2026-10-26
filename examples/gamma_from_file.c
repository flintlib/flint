#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include <flint/flint.h>
#include <flint/acb.h>

const slong bin_prec = 128; // binary precision for computation
const long dec_prec = 38; // decimal digits

void acb_set_str(acb_t z, const char *real_str, const char *imag_str) {

    arb_t real_part, imag_part;
    arb_init(real_part);
    arb_init(imag_part);

    // Set the real and imaginary parts from strings
    if (arb_set_str(real_part, real_str, 128) == 1) {
        fprintf(stderr, "Error setting real part from string: %s\n", real_str);
        arb_clear(real_part);
        arb_clear(imag_part);
        return;
    }

    if (arb_set_str(imag_part, imag_str, 128) == 1) {
        fprintf(stderr, "Error setting imaginary part from string: %s\n", imag_str);
        arb_clear(real_part);
        arb_clear(imag_part);
        return;
    }

    // Combine into acb_t
    acb_set_arb_arb(z, real_part, imag_part);
    arb_clear(real_part);
    arb_clear(imag_part);
}

void parse_complex_number_str(const char *input_line, char *real_part, char *imag_part) {
    // Find the first space character
    char *space_pos = strchr(input_line, ' ');
    if (space_pos) {
        // Get the real part
        strncpy(real_part, input_line, space_pos - input_line);
        real_part[space_pos - input_line] = '\0'; // Null-terminate

        // Get the imaginary part
        // Skip the space and copy the rest
        strcpy(imag_part, space_pos + 1);

        // Remove spaces from the imaginary part
        char *src = imag_part, *dst = imag_part;
        while (*src) {
            if (*src != ' ') {
                *dst++ = *src; // Copy non-space characters
            }
            src++;
        }
        *dst = '\0'; // Null terminate the modified imaginary part

        // Remove all spaces and the '*I' from the imaginary part
        char *imaginary_end = strstr(imag_part, "*I");
        if (imaginary_end) {
            *imaginary_end = '\0'; // Cut off the '*I'
        }

    } else {
        fprintf(stderr, "Error: No space found in input line: '%s'\n", input_line);
    }
}

void compute_gamma_from_file(const char *filename) {
    FILE *file = fopen(filename, "r");
    if (file == NULL) {
        perror("Error opening file");
        return;
    }

    char input_line[2 * dec_prec + 20];
    char real_str[dec_prec + 10];
    char imag_str[dec_prec + 10];

    acb_t z, result;
    acb_init(z);
    acb_init(result);
    unsigned long records = 0;

    while (fgets(input_line, sizeof(input_line), file)) {

        memset(real_str, 0, sizeof(real_str));
        memset(imag_str, 0, sizeof(imag_str));

        // Find and remove newline character
        char *newline_pos = strchr(input_line, '\n');
        if (newline_pos) {
            *newline_pos = '\0'; // Replace newline with null terminator
        }

        // Print the raw input line for debugging
/*
        printf("Raw Input Line: '%s'\n", input_line);
*/
        if (strlen(input_line)) {
            parse_complex_number_str(input_line, real_str, imag_str);
        }

        // Print parsed values for debugging
/*
        printf("Real:'%s'\n", real_str);
        printf("Imag:'%s'\n", imag_str);
*/
        acb_set_str(z, real_str, imag_str);
        acb_gamma(result, z, bin_prec);  // Compute Gamma(z)

        // Print result
        // Check for a minus sign in the imaginary part
        printf("Gamma(%s%s*I) = ", real_str, imag_str);
        acb_printd(result, 38); // Print with 38 decimal places
        printf("\n");

        records++; // Increment records for successfully processed line

    }

    printf("Read: %lu records\n", records);

    acb_clear(z);
    acb_clear(result);
    fclose(file);
}

int main(int argc, char *argv[]) {
    if (argc != 2) {
        fprintf(stderr, "Usage: %s <filename>\n", argv[0]);
        return EXIT_FAILURE;
    }

    clock_t start_time = clock(); // Start timing

    compute_gamma_from_file(argv[1]);

    clock_t end_time = clock(); // End timing
    double time_taken = (double)(end_time - start_time) / CLOCKS_PER_SEC; // Convert to seconds
    printf("Time taken: %.6f seconds\n", time_taken);

    return EXIT_SUCCESS;
}

/* Input file format: <real part 38 digits> <sign of imag> <abs(imag part 38 digits)>*I
------------------------------------------------------------------------------------
-7.6863655401393771171569824218750000000 - 3.3994659781455993652343750000000000000*I
-5.6718241423368453979492187500000000000 - 2.4417644459754228591918945312500000000*I
etc.
------------------------------------------------------------------------------------
*/

