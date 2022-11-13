import java.lang.Math;

/**** on va ici implémenter la transformée de Fourier rapide 1D ****/

public class FFT_1D_part2 {

	// "combine" c1 et c2 selon la formule vue en TD
	// c1 et c2 sont de même taille
	// la taille du résultat est le double de la taille de c1
	public static CpxTab combine(CpxTab c1, CpxTab c2) {
		assert (c1.taille() == c2.taille())
				: "combine: c1 et c2 ne sont pas de même taille, taille c1=" + c1.taille() + " taille c2=" + c2.taille();
		// A FAIRE
		int n = 2 * c1.taille(); // On stock la taille de c1
		/*
		 * CpxTab c va stocker la combinaison de c1 et c2
		 * c[k] = c1[k] + w[k] et c[k+n/2] = c1[k+n/2] + w[k+n/2]
		 */
		CpxTab c = new CpxTab(n);
		// La boucle for suivante permet de remplir c
		for (int k = 0; k < n / 2; k++) {
			// on calcul les partie réelles et imaginaire de c
			double c_k_reel = c1.get_p_reel(k) + Math.cos((2 * k * Math.PI) / n) * c2.get_p_reel(k)
					- Math.sin((2 * k * Math.PI) / n) * c2.get_p_imag(k);
			double c_k_imag = c1.get_p_imag(k) + Math.sin((2 * k * Math.PI) / n) * c2.get_p_reel(k)
					+ Math.cos((2 * k * Math.PI) / n) * c2.get_p_imag(k);
			double c_kn2_reel = c1.get_p_reel(k) - Math.cos((2 * k * Math.PI) / n) * c2.get_p_reel(k)
					+ Math.sin((2 * k * Math.PI) / n) * c2.get_p_imag(k);
			double c_kn2_imag = c1.get_p_imag(k) - Math.cos((2 * k * Math.PI) / n) * c2.get_p_imag(k)
					- Math.sin((2 * k * Math.PI) / n) * c2.get_p_reel(k);
			// on stosk ces partie réelles et imaginaire dans c
			c.set_p_reel(k, c_k_reel);
			c.set_p_imag(k, c_k_imag);
			c.set_p_reel(k + n / 2, c_kn2_reel);
			c.set_p_imag(k + n / 2, c_kn2_imag);
		}
		return c;
	}

	// renvoie la TFD d'un tableau de complexes
	// la taille de x doit être une puissance de 2
	public static CpxTab FFT(CpxTab x) {
		// A FAIRE : Test d'arrêt
		if (x.taille() == 1)
			return x;

		assert (x.taille() % 2 == 0) : "FFT: la taille de x doit être une puissance de 2";

		// A FAIRE : Décomposition en "pair" et "impair" et appel récursif
		int n = x.taille() / 2;
		CpxTab x_pair = new CpxTab(n);
		CpxTab x_impair = new CpxTab(n);
		int i = 0;
		int j = 0;
		while (i < n || j < n) {
			for (int k = 0; k < 2 * n; k++) {
				if (k % 2 == 0) {
					x_pair.set_p_reel(i, x.get_p_reel(k));
					x_pair.set_p_imag(i, x.get_p_imag(k));
					i++;
				} else {
					x_impair.set_p_reel(j, x.get_p_reel(k));
					x_impair.set_p_imag(j, x.get_p_imag(k));
					j++;
				}
			}
		}
		CpxTab y_pair = FFT(x_pair);
		CpxTab y_impair = FFT(x_impair);
		CpxTab y = combine(y_pair, y_impair);
		return y;
	}

	// renvoie la TFD d'un tableau de réels
	// la taille de x doit être une puissance de 2
	public static CpxTab FFT(double[] x) {
		return FFT(new CpxTab(x));
	}

	// renvoie la transformée de Fourier inverse de y
	public static CpxTab FFT_inverse(CpxTab y) {
		// A FAIRE
		int n = y.taille();
		CpxTab c1 = (FFT(y.conjugue())).conjugue();
		CpxTab y_FFT_inverse = new CpxTab(n);
		for (int i = 0; i < n; i++) {
			y_FFT_inverse.set_p_reel(i, c1.get_p_reel(i) / n);
			y_FFT_inverse.set_p_imag(i, c1.get_p_imag(i) / n);
		}
		return y_FFT_inverse;
	}

	// calcule le produit de deux polynômes en utilisant la FFT
	// tab1 et tab2, sont les coefficients de ces polynômes
	// CpxTab sera le tableau des coefficients du polynôme produit (purement réel)
	public static CpxTab multiplication_polynome_viaFFT(double[] tab1, double[] tab2) {

		// on commence par doubler la taille des vecteurs en rajoutant des zéros à la
		// fin (cf TD)
		double[] t1 = new double[2 * tab1.length], t2 = new double[2 * tab1.length];
		for (int i = 0; i < tab1.length; i++) {
			t1[i] = tab1[i];
			t2[i] = tab2[i];
		}

		// A COMPLETER !!

		return null;
	}

	// renvoie un tableau de réels aléatoires
	// utile pour tester la multiplication de polynômes
	public static double[] random(int n) {
		double[] t = new double[n];

		for (int i = 0; i < n; i++)
			t[i] = Math.random();
		return t;
	}

	// effectue la multiplication de polynômes représentés par coefficients
	// p1, p2 les coefficients des deux polynômes P1 et P2
	// renvoie les coefficients du polynôme P1*P2
	private static double[] multiplication_polynome_viaCoeff(double[] p1, double[] p2) {

		int n = p1.length + p2.length - 1;
		double a, b;
		double[] out = new double[n];
		for (int k = 0; k < n; k++) {
			for (int i = 0; i <= k; i++) {
				a = (i < p1.length) ? p1[i] : 0;
				b = (k - i < p2.length) ? p2[k - i] : 0;
				out[k] += a * b;
			}
		}
		return out;
	}

	// affiche un tableau de réels
	private static void afficher(double[] t) {
		System.out.print("[");
		for (int k = 0; k < t.length; k++) {
			System.out.print(t[k]);
			if (k < (t.length - 1))
				System.out.print(" ");
		}
		System.out.println("]");
	}

	public static void main(String[] args) {
		/* PARTI 2 */
		/* Exo 1: Tests */
			int n = 16;
		
			/* TFD d'un vecteur constant a = (a, ..., a) */
			System.out.println("-----------------------------------------------------");
			System.out.println("TFD d'un vecteur constant a = (a, ..., a)");
			double []k = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12};
			for (double a : k){
				double[] t6 = new double[n];
				for (int i = 0; i < n; i++)
					t6[i] = a;
				System.out.println("	Pour a="+a+" : "+FFT(new CpxTab(t6)).toString());
			}
			/* TFD d'une sinusoïde pure */
			System.out.println("-----------------------------------------------------");
			System.out.println("TFD d'une sinusoïde pure a = (cos((2πk/n)*0), cos((2πk/n)*1), ..., cos((2πk/n)*("+(n-1)+")))");
			double[] t7 = new double[n];
			for (double a : k){
				for (int i = 0; i < n; i++) {
					t7[i] = Math.cos(2 * Math.PI * a * i / n);
				}
				System.out.println("	Pour k="+a+" : "+FFT(new CpxTab(t7)).toString());
			}
		
			/* TFD somme de 2 sinusoïdes pures */
			System.out.println("-----------------------------------------------------");
			System.out.println("TFD somme de 2 sinusoïdes pures");
			for (int i = 0; i < n; i++) {
				t7[i] = Math.cos(2 * Math.PI * i / n) + 0.5 * Math.cos(2 * Math.PI * 3 * i / n);
			}
			System.out.println("	"+(FFT(new CpxTab(t7)).toString()));
		
			/* TFD somme de 2 sinusoïdes pures et d'une constante */
			System.out.println("-----------------------------------------------------");
			System.out.println("TFD somme de 2 sinusoïdes pures et d'une constante");
			for (int i = 0; i < n; i++) {
				t7[i] = 4 + 2 * Math.sin(2 * Math.PI * 2 * i / n) + 0.5 * Math.cos(2 * Math.PI * 7 * i / n);
			}
			System.out.println("	"+FFT(new CpxTab(t7)).toString());
	}
}