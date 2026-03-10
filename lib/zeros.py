import flint

flint.ctx.prec = 96

def get_zero(n):

    pi = flint.arb.pi()

    lamb = ((n-0.125)*pi).sqrt()

    mu = (2*(2*pi).sqrt()*lamb).log()

    term2 = mu * 0.25 * lamb**-1

    term3 = (1/16) * (1 - mu + 0.5 * mu**2) * (lamb**-3)

    xn = - lamb + term2 - term3
    yn = lamb + term2 + term3

    return flint.acb(xn,yn)


def newton(z_guess, it=15):
    z=z_guess
    deriv_const = -2 / flint.arb.pi().sqrt()

    for _ in range(it):
        #w=(-z**2).exp() * (-flint.acb(0,1)*z).erfc()
        #w2=-2*z*w+2*flint.acb(0,1) / flint.acb.pi().sqrt()
        #step = w/w2 

        f=z.erfc()
        f2=deriv_const*(-z**2).exp()
        step=f/f2

        if step.abs_upper() < (z.abs_lower()*flint.arb("1e-24")):
            break

        z=z-step

    return z



def main():
    
    for n in range(1, 11):
        z_approx = get_zero(n)
        z_exact = newton(z_approx)
        print(f"{n:<} | {z_exact.str(24)}")

if __name__ == "__main__":
    main()
