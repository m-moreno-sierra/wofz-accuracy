#include <stdio.h>
#include <flint/acb.h>
#include <flint/acb_hypgeom.h>

void find_y(arb_t result, arb_t x, arb_t C, arb_t tol) {
    slong prec = 64;  // Set precision
    arb_t y, y_low, y_high, mid, real_part, diff, abs_diff;
    acb_t z, w;

    arb_init(y);
    arb_init(y_low);
    arb_init(y_high);
    arb_init(mid);
    arb_init(real_part);
    arb_init(diff);
    arb_init(abs_diff);

    acb_init(z);
    acb_init(w);

    // Initial bounds for y
    arb_set_si(y_low, -10);  // Lower bound
    arb_set_si(y_high, 10);  // Upper bound

    while (1) {
        // Midpoint
        arb_add(mid, y_low, y_high, prec);
        arb_mul_2exp_si(mid, mid, -1);

        // Evaluate z = x + i * mid
        acb_set_arb_arb(z, x, mid);
        acb_hypgeom_erfc(w, z, prec);  // Compute w(z)
        acb_mul(w, w, z, prec);        // Faddeeva function: w(z) = e^(-z^2) * erfc(-i * z)

        // Get real part
        acb_get_real(real_part, w);

        // Compute diff = Re(w(z)) - C
        arb_sub(diff, real_part, C, prec);

        // Compute absolute value and check tolerance
        arb_abs(abs_diff, diff);  // abs_diff = |diff|
        if (arb_lt(abs_diff, tol)) {
            arb_set(result, mid);  // y found
            break;
        }

        // Update bounds
        if (arb_is_positive(diff)) {
            arb_set(y_high, mid);  // diff > 0, y is smaller
        } else {
            arb_set(y_low, mid);   // diff < 0, y is larger
        }
    }

    // Clear memory
    arb_clear(y);
    arb_clear(y_low);
    arb_clear(y_high);
    arb_clear(mid);
    arb_clear(real_part);
    arb_clear(diff);
    arb_clear(abs_diff);
    acb_clear(z);
    acb_clear(w);
}

int main() {
    arb_t x, C, tol, y;
    arb_init(x);
    arb_init(C);
    arb_init(tol);
    arb_init(y);

    // Example values
    arb_set_d(x, 1.0);
    arb_set_d(C, 0);
    arb_set_d(tol, 1e-10);

    find_y(y, x, C, tol);

    // Print result
    printf("y = ");
    arb_printn(y, 10, 0);
    printf("\n");

    // Clear memory
    arb_clear(x);
    arb_clear(C);
    arb_clear(tol);
    arb_clear(y);

    return 0;
}
