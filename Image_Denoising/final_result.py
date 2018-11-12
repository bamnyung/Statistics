import os

import matplotlib.pyplot as plt
import numpy as np


class Denoising(object):
    def __init__(self, img, noise=True):
        self.img = img
        self.noise = noise
        self.is_gray = False
        self.is_jpg = False
        self.method = None

    @staticmethod
    def px(u):
        return np.roll(u, -1, axis=1) - u

    @staticmethod
    def px2(u):
        return np.roll(u, -1, axis=1) - np.roll(u, 1, axis=1)

    @staticmethod
    def py(u):
        return np.roll(u, -1, axis=0) - u

    @staticmethod
    def py2(u):
        return np.roll(u, -1, axis=0) - np.roll(u, 1, axis=0)

    @staticmethod
    def mx(u):
        return u - np.roll(u, 1, axis=1)

    @staticmethod
    def my(u):
        return u - np.roll(u, 1, axis=0)

    @staticmethod
    def minmod(a, b):
        return (np.sign(a) + np.sign(b)) / 2 * np.minimum(abs(a), abs(b))

    @staticmethod
    def H(mtx1, mtx2):
        # p=1
        return 1 / (mtx1 ** 2 + mtx2 ** 2 + 1.0E-06 ** 2) ** 0.5

    def rgb2gray(self):
        rgb = self.img
        self.is_gray = True
        if rgb.ndim == 2:
            self.img = rgb
        else:
            self.img = np.dot(rgb[:, :, :3], [0.299, 0.587, 0.114])

    def denoise_level(self, tol=1.0E-06, max_iter_num=200, sigma=0.1, dt=0.01):
        self.method = 'level'
        self.rgb2gray()
        self.adjust_image()
        u = self.img
        n, m = u.shape[:2]
        h1 = 1 / n
        h2 = 1 / m

        i = 0
        while i < max_iter_num:
            u_old = u

            dx = np.divide(
                self.px(u),
                np.sqrt(np.square(self.px(u)) +
                        np.square(self.minmod(self.py(u), self.my(u)))) +
                tol ** 2
            )
            dy = np.divide(
                self.py(u),
                np.sqrt(np.square(self.py(u)) +
                        np.square(self.minmod(self.px(u), self.mx(u)))) +
                tol ** 2
            )

            normu = np.sqrt(np.square(self.px(u)) + np.square(self.py(u)))
            lam = -np.sqrt(h1 ** 2 + h2 ** 2) / (2 * sigma ** 2) * (sum(sum(
                normu -
                np.divide(np.multiply(self.px(self.img),
                                      self.px(u)), normu + tol ** 2) -
                np.divide(np.multiply(self.py(self.img),
                                      self.py(u)), normu + tol ** 2)
            )))

            u = u + dt/h1*self.mx(dx) + dt/h2*self.my(dy) - dt*lam*(u-self.img)

            error = np.linalg.norm(u - u_old) / np.sqrt(n * m)

            if i == 0:
                err_init = error
                err_prev = error
            else:
                if np.abs(err_prev - error) < tol * err_init:
                    break
                else:
                    err_prev = error

            i += 1
        self.denoised = u

    def denoise_texture(self, tol=1.0E-06, max_iter_num=200, lam=0.1, mu=0.1):
        self.method = 'texture'
        self.rgb2gray()
        self.adjust_image()
        u = self.img
        n, m = u.shape[:2]
        h1 = 1 / n
        h2 = 1 / m
        fx = self.px2(u)
        fy = self.py2(u)
        g1 = -fx / (2 * lam * (fx ** 2 + fy ** 2) ** 0.5 + tol ** 2)
        g2 = -fy / (2 * lam * (fx ** 2 + fy ** 2) ** 0.5 + tol ** 2)

        i = 0

        while i < max_iter_num:
            u_old = u

            c1 = 1 / ((self.px(u) / h1) ** 2 + (
                self.py2(u) / (2 * h2)) ** 2 + tol ** 2) ** 0.5
            c2 = 1 / ((self.mx(u) / h1) ** 2 +
                      (self.py2(np.roll(u, 1, axis=1)) / (2 * h2)) ** 2 +
                      tol ** 2) ** 0.5
            c3 = 1 / ((self.px2(u) / h1) ** 2 + (self.py(u) / h2) ** 2 +
                      tol ** 2) ** 0.5
            c4 = 1 / ((self.px2(np.roll(u, 1, axis=0)) / h1) ** 2 +
                      (self.my(u) / h2) ** 2 + tol ** 2) ** 0.5

            u = np.multiply(
                1 / (1 + (c1 + c2 + c3 + c4) / (2 * lam * h1 * h2)),
                self.img - self.px2(g1) / (2 * h1) - self.py2(g2) / (2 * h2) +
                (np.multiply(c1, np.roll(u, -1, axis=1)) +
                 np.multiply(c2, np.roll(u, 1, axis=1)) +
                 np.multiply(c3, np.roll(u, -1, axis=0)) +
                 np.multiply(c4, np.roll(u, 1, axis=0))) / (2 * lam * h1 * h2)
            )

            g1 = np.multiply(
                2 * lam / (mu * self.H(g1, g2) + 4 * lam / (h1 ** 2)),
                self.px2(u) / (2 * h1) - self.px2(self.img) / (2 * h1) +
                self.px2(g1) / (h1 ** 2) +
                (self.px(np.roll(g2, -1, axis=0)) - self.px(g2)) / (2*h1**2) +
                (self.my(g2) - self.my(np.roll(g2, 1, axis=1))) / (2*h2**2)
            )

            g2 = np.multiply(
                2 * lam / (mu * self.H(g1, g2) + 4 * lam / (h2 ** 2)),
                self.py2(u) / (2 * h2) - self.py2(self.img) / (2 * h2) +
                self.py2(g2) / (h2 ** 2) +
                (self.px(np.roll(g1, -1, axis=0)) - self.px(g1)) / (2*h1**2) +
                (self.my(g1) - self.my(np.roll(g1, 1, axis=1))) / (2*h2**2)
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

            i += 1
        self.denoised = u

    @staticmethod
    def _denoise_tv_chambolle_nd(im, weight=0.1, eps=2.e-4, n_iter_max=200):
        ndim = im.ndim
        p = np.zeros((im.ndim,) + im.shape, dtype=im.dtype)
        g = np.zeros_like(p)
        d = np.zeros_like(im)
        i = 0
        while i < n_iter_max:
            if i > 0:
                # d will be the (negative) divergence of p
                d = -p.sum(0)
                slices_d = [slice(None), ] * ndim
                slices_p = [slice(None), ] * (ndim + 1)
                for ax in range(ndim):
                    slices_d[ax] = slice(1, None)
                    slices_p[ax + 1] = slice(0, -1)
                    slices_p[0] = ax
                    d[slices_d] += p[slices_p]
                    slices_d[ax] = slice(None)
                    slices_p[ax + 1] = slice(None)
                out = im + d
            else:
                out = im
            E = (d ** 2).sum()

            # g stores the gradients of out along each axis
            # e.g. g[0] is the first order finite difference along axis 0
            slices_g = [slice(None), ] * (ndim + 1)
            for ax in range(ndim):
                slices_g[ax + 1] = slice(0, -1)
                slices_g[0] = ax
                g[slices_g] = np.diff(out, axis=ax)
                slices_g[ax + 1] = slice(None)

            norm = np.sqrt((g ** 2).sum(axis=0))[np.newaxis, ...]
            E += weight * norm.sum()
            tau = 1. / (2. * ndim)
            norm *= tau / weight
            norm += 1.
            p -= tau * g
            p /= norm
            E /= float(im.size)
            if i == 0:
                E_init = E
                E_previous = E
            else:
                if np.abs(E_previous - E) < eps * E_init:
                    break
                else:
                    E_previous = E
            i += 1
        return out

    def denoise_tv_chambolle(self, weight=0.1, eps=2.e-4, n_iter_max=200):
        self.method = 'chambolle'
        self.adjust_image()
        im = self.img

        if self.is_gray:
            out = self._denoise_tv_chambolle_nd(im, weight, eps, n_iter_max)
        else:
            out = np.zeros_like(im)
            for c in range(im.shape[-1]):
                out[..., c] = self._denoise_tv_chambolle_nd(im[..., c],
                                                            weight, eps,
                                                            n_iter_max)
        self.denoised = out

    def adjust_image(self):
        img = self.img
        if img.ndim != 3:
            self.is_gray = True
            if sum(sum(img <= 1)) != img.size:
                self.is_jpg = True
        else:
            if img[:, :, 0].size not in [sum(sum(img[:, :, i] <= 1)) for i in
                                     range(
                    img.shape[2])]:
                self.is_jpg = True
        if self.is_jpg:
            img = img / 255
            self.img = img
        if not self.noise:
            img = img + 0.1 * np.random.standard_normal(self.img.shape)
            img = np.clip(img, 0, 1)
            self.img = img

    def plot_image(self):
        if self.is_gray:
            plt.subplot(121)
            plt.imshow(self.img, cmap='gray')
            plt.title('real image')

            plt.subplot(122)
            plt.imshow(self.denoised, cmap='gray')
            plt.title('changed')
        else:
            plt.subplot(121)
            plt.imshow(self.img)
            plt.title('real image')

            plt.subplot(122)
            plt.imshow(self.denoised)
            plt.title('changed')

        if not os.path.exists('result'):
            os.makedirs('result')

        if self.is_gray:
            plt.imsave('result/image_gray.jpg', self.img, cmap='gray')
        else:
            plt.imsave('result/image_real.jpg', self.img)
        plt.imsave('result/changed_{}.jpg'.format(self.method), self.denoised,
                   cmap='gray')
        plt.suptitle('image denoising with {} method'.format(self.method))
        plt.show()

if __name__ == "__main__":
    fName = ""
    while True:
        fName = input("Please enter the name of your file: ")
        if os.access(fName, os.R_OK) : break
        else :
            print("Your file does not exist in the current directory.")

    image = plt.imread(fName)

    img1 = Denoising(image, noise=False)
    img1.denoise_level(max_iter_num=10, dt=0.00001)
    img1.plot_image()
    img2 = Denoising(image, noise=False)
    img2.denoise_texture(tol=1.0E-06, max_iter_num=3, lam=0.1, mu=0.01)
    img2.plot_image()
    img3 = Denoising(image, noise=False)
    img3.denoise_tv_chambolle()
    img3.plot_image()