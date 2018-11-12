import matplotlib.pyplot as plt
import numpy as np

import matplotlib.pyplot as plt
import numpy as np


def denoise_level(img, tol=1.0E-06, max_iter_num=200, sigma=30, dt=0.001):
    u = img
    n, m = u.shape[:2]
    h1 = 1/n
    h2 = 1/m

    i = 0
    while i < max_iter_num:
        u_old = u

        dx = np.divide(
            px(u),
            np.sqrt(np.square(px(u)) + np.square(minmod(py(u), my(u)))) +
            tol**2
        )
        dy = np.divide(
            py(u),
            np.sqrt(np.square(py(u)) + np.square(minmod(px(u), mx(u)))) +
            tol**2
        )

        normu = np.sqrt(np.square(px(u)) + np.square(py(u)))
        lam = -h1/(2 * sigma ** 2) * (sum(sum(
            normu -
            np.divide(np.multiply(px(img), px(u)), normu + tol**2) -
            np.divide(np.multiply(py(img), py(u)), normu + tol**2)
        )))

        u = u + dt/h1*mx(dx) + dt/h2*my(dy) - dt*lam*(u - img)

        error = np.linalg.norm(u - u_old) / np.sqrt(n*m)

        if i == 0:
            err_init = error
            err_prev = error
        else:
            if np.abs(err_prev - error) < tol * err_init:
                break
            else:
                err_prev = error

        # don't forget to update iterator
        i += 1
        print(i, error)

    return u, i


def denoise_texture(img, tol=1.0E-06, max_iter_num=200, lam=0.1, mu=0.1):
    u = img
    n, m = u.shape[:2]
    h1 = 1 / n
    h2 = 1 / m
    fx = px2(u)
    fy = py2(u)
    g1 = -fx/(2*lam*(fx**2 + fy**2)**0.5 + tol**2)
    g2 = -fy/(2*lam*(fx**2 + fy**2)**0.5 + tol**2)

    i = 0

    while i < max_iter_num:
        u_old = u

        c1 = 1 / ((px(u)/h1)**2 + (py2(u)/(2*h2))**2 + tol**2)**0.5
        c2 = 1 / ((mx(u)/h1)**2 +
                  (py2(np.roll(u, 1, axis=1))/(2*h2))**2 +
                  tol**2)**0.5
        c3 = 1 / ((px2(u)/h1)**2 + (py(u)/h2)**2 + tol**2)**0.5
        c4 = 1 / ((px2(np.roll(u, 1, axis=0))/h1)**2 +
                  (my(u)/h2)**2 + tol**2)**0.5

        u = np.multiply(
            1/(1+(c1+c2+c3+c4)/(2*lam*h1*h2)),
            img - px2(g1)/(2*h1) - py2(g2)/(2*h2) +
            (np.multiply(c1, np.roll(u, -1, axis=1)) +
             np.multiply(c2, np.roll(u, 1, axis=1)) +
             np.multiply(c3, np.roll(u, -1, axis=0)) +
             np.multiply(c4, np.roll(u, 1, axis=0))) / (2*lam*h1*h2)
        )

        g1 = np.multiply(
            2*lam/(mu*H(g1, g2) + 4*lam/(h1**2)),
            px2(u)/(2*h1) - px2(img)/(2*h1) + px2(g1)/(h1**2) +
            (px(np.roll(g2, -1, axis=0)) - px(g2))/(2*h1**2) +
            (my(g2) - my(np.roll(g2, 1, axis=1)))/(2*h2**2)
        )

        g2 = np.multiply(
            2*lam/(mu*H(g1, g2) + 4*lam/(h2**2)),
            py2(u)/(2*h2) - py2(img)/(2*h2) + py2(g2)/(h2**2) +
            (px(np.roll(g1, -1, axis=0)) - px(g1))/(2*h1**2) +
            (my(g1) - my(np.roll(g1, 1, axis=1)))/(2*h2**2)
        )

        error = np.linalg.norm(u - u_old) / np.sqrt(n * m)

        if i == 0:
            err_init = error
            err_prev = error
        else:
            if np.abs(err_prev - error) < tol * err_init:
                break
            else:
                err_prev = error

        # don't forget to update iterator
        i += 1
        print(i, error)

    return u, i


def px(u):
    return np.roll(u, -1, axis=1) - u


def px2(u):
    return np.roll(u, -1, axis=1) - np.roll(u,1,axis=1)


def py(u):
    return np.roll(u, -1, axis=0) - u


def py2(u):
    return np.roll(u, -1, axis=0) - np.roll(u,1,axis=0)


def mx(u):
    return u - np.roll(u, 1, axis=1)


def my(u):
    return u - np.roll(u, 1, axis=0)


def minmod(a, b):
    return (np.sign(a) + np.sign(b))/2 * np.minimum(abs(a), abs(b))


def H(mtx1, mtx2):
    # p=1
    return 1 / (mtx1**2+mtx2**2 + 1.0E-06**2)**0.5


def rgb2gray(rgb):
    if rgb.ndim == 2:
        return rgb
    else:
        return np.dot(rgb[:, :, :3], [0.299, 0.587, 0.114])

if __name__ == "__main__":
    # img_real = ascent()
    img_real = plt.imread('image.png')
    img_gray = rgb2gray(img_real)

    if sum(sum(img_gray <= 1)) == img_gray.size:
        img_gray *= 255
    sigma = 30
    img_noise = img_gray + sigma * np.random.standard_normal(img_gray.shape)
    img_noise = np.clip(img_noise, 0, 255)

    # img_denoise, i = denoise_level(img_noise, max_iter_num=10, dt=0.001)
    img_denoise, i = denoise_texture(
        img_noise, tol=1.0E-06, max_iter_num=3, lam=0.1, mu=0.01)

    plt.subplot(141)
    plt.imshow(img_real)
    plt.title('real image')

    plt.subplot(142)
    plt.imshow(img_gray, cmap="gray")
    plt.title('gray scale')

    plt.subplot(143)
    plt.imshow(img_noise, cmap="gray")
    plt.title('noisy')

    plt.subplot(144)
    plt.imshow(img_denoise, cmap="gray")
    plt.title('denoising')

    plt.show()
