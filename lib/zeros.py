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

def main():
    
    for n in range(1, 6):
        z_approx = get_zero(n)
        print(f"{n:<3} | {z_approx.str(20)}")

if __name__ == "__main__":
    main()
