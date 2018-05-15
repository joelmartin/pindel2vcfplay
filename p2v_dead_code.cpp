// ################################
// log\binom{n}{k}
static double lbinom(int n, int k) {
	if (k == 0 || n == k) {
		return 0;
	}
	return lgamma(n+1) - lgamma(k+1) - lgamma(n-k+1);
}

// n11  n12  | n1_
// n21  n22  | n2_
//-----------+----
// n_1  n_2  | n

// hypergeometric distribution
static double hypergeo(int n11, int n1_, int n_1, int n) {
	return exp(lbinom(n1_, n11) + lbinom(n-n1_, n_1-n11) - lbinom(n, n_1));
}

typedef struct {
	int n11, n1_, n_1, n;
	double p;
} hgacc_t;

// incremental version of hypergenometric distribution
static double hypergeo_acc(int n11, int n1_, int n_1, int n, hgacc_t *aux) {
	if (n1_ || n_1 || n) {
		aux->n11 = n11;
		aux->n1_ = n1_;
		aux->n_1 = n_1;
		aux->n = n;
	}
	else {   // then only n11 changed; the rest fixed
		if (n11%11 && n11 + aux->n - aux->n1_ - aux->n_1) {
			if (n11 == aux->n11 + 1) { // incremental
				aux->p *= (double)(aux->n1_ - aux->n11) / n11
				          * (aux->n_1 - aux->n11) / (n11 + aux->n - aux->n1_ - aux->n_1);
				aux->n11 = n11;
				return aux->p;
			}
			if (n11 == aux->n11 - 1) { // incremental
				aux->p *= (double)aux->n11 / (aux->n1_ - n11)
				          * (aux->n11 + aux->n - aux->n1_ - aux->n_1) / (aux->n_1 - n11);
				aux->n11 = n11;
				return aux->p;
			}
		}
		aux->n11 = n11;
	}
	aux->p = hypergeo(aux->n11, aux->n1_, aux->n_1, aux->n);
	return aux->p;
}

double kt_fisher_exact(int n11, int n12, int n21, int n22, double *_left, double *_right, double *two) {
	int i, j, max, min;
	double p, q, left, right;
	hgacc_t aux;
	int n1_, n_1, n;

	n1_ = n11 + n12;
	n_1 = n11 + n21;
	n = n11 + n12 + n21 + n22; // calculate n1_, n_1 and n
	max = (n_1 < n1_) ? n_1 : n1_; // max n11, for right tail
	min = n1_ + n_1 - n;
	if (min < 0) {
		min = 0;   // min n11, for left tail
	}
	*two = *_left = *_right = 1.;
	if (min == max) {
		return 1.;   // no need to do test
	}
	q = hypergeo_acc(n11, n1_, n_1, n, &aux); // the probability of the current table
	// left tail
	p = hypergeo_acc(min, 0, 0, 0, &aux);
	for (left = 0., i = min + 1; p < 0.99999999 * q; ++i) { // loop until underflow
		left += p, p = hypergeo_acc(i, 0, 0, 0, &aux);
	}
	--i;
	if (p < 1.00000001 * q) {
		left += p;
	}
	else {
		--i;
	}
	// right tail
	p = hypergeo_acc(max, 0, 0, 0, &aux);
	for (right = 0., j = max - 1; p < 0.99999999 * q; --j) { // loop until underflow
		right += p, p = hypergeo_acc(j, 0, 0, 0, &aux);
	}
	++j;
	if (p < 1.00000001 * q) {
		right += p;
	}
	else {
		++j;
	}
	// two-tail
	*two = left + right;
	if (*two > 1.) {
		*two = 1.;
	}
	// adjust left and right
	if (abs(i - n11) < abs(j - n11)) {
		right = 1. - left + q;
	}
	else {
		left = 1.0 - right + q;
	}
	*_left = left;
	*_right = right;
	return q;
}

#define EXACT_TEST_BIAS 0.00000000000000000000000010339757656912845935892608650874535669572651386260986328125

double fisher22_1sided(uint32_t m11, uint32_t m12, uint32_t m21, uint32_t m22, uint32_t m11_is_greater_alt) {
	double cur_prob = EXACT_TEST_BIAS;
	double left_prob = cur_prob;
	double right_prob = 0;
	uint32_t uii;
	double cur11;
	double cur12;
	double cur21;
	double cur22;
	double preaddp;
	// Ensure m11 <= m22 and m12 <= m21.
	if (m12 > m21) {
		uii = m12;
		m12 = m21;
		m21 = uii;
	}
	if (m11 > m22) {
		uii = m11;
		m11 = m22;
		m22 = uii;
	}
	// Flipping m11<->m12 and m21<->m22 also flips the direction of the
	// alternative hypothesis.  So we flip on m11-is-greater alternative
	// hypothesis here to allow the rest of the code to assume m11-is-less.
	if (m11_is_greater_alt) {
		uii = m11;
		m11 = m12;
		m12 = uii;
		uii = m21;
		m21 = m22;
		m22 = uii;
	}
	cur11 = m11;
	cur12 = m12;
	cur21 = m21;
	cur22 = m22;
	if ((((uint64_t)m11) * m22) >= (((uint64_t)m12) * m21)) {
		// starting right of (or at) center, p > 0.5
		// 1. left_prob = sum leftward to precision limit
		// 2. total_prob := left_prob
		// 3. total_prob += sum rightward to total_prob precision limit
		// return left_prob / total_prob
		while (cur11 > 0.5) {
			cur12 += 1;
			cur21 += 1;
			cur_prob *= (cur11 * cur22) / (cur12 * cur21);
			cur11 -= 1;
			cur22 -= 1;
			preaddp = left_prob;
			left_prob += cur_prob;
			if (left_prob <= preaddp) {
				break;
			}
			if (left_prob >= 1.0) {
				// Probability mass of our starting table was represented as 2^{-83},
				// so this would mean the left probability mass partial sum is greater
				// than 2^83 times that.  In which case the final p-value will
				// be indistinguishable from 1 at 53-bit precision if our input just
				// had 32-bit integers.  (Yes, the constant can be reduced.)
				return 1;
			}
		}
		cur11 = m11;
		cur12 = m12;
		cur21 = m21;
		cur22 = m22;
		cur_prob = EXACT_TEST_BIAS;
		right_prob = left_prob; // actually total_prob
		while (cur12 > 0.5) {
			cur11 += 1;
			cur22 += 1;
			cur_prob *= (cur12 * cur21) / (cur11 * cur22);
			cur12 -= 1;
			cur21 -= 1;
			preaddp = right_prob;
			right_prob += cur_prob;
			if (right_prob <= preaddp) {
				break;
			}
		}
		return left_prob / right_prob;
	}
	else {
		// starting left of center, p could be small
		// 1. right_prob = sum rightward to precision limit
		// 2. left_prob = sum leftward to left_prob precision limit
		// return left_prob / (left_prob + right_prob)
		while (cur12 > 0.5) {
			cur11 += 1;
			cur22 += 1;
			cur_prob *= (cur12 * cur21) / (cur11 * cur22);
			cur12 -= 1;
			cur21 -= 1;
			preaddp = right_prob;
			right_prob += cur_prob;
			if (right_prob == INFINITY) {
				return 0;
			}
			if (right_prob <= preaddp) {
				break;
			}
		}
		cur11 = m11;
		cur12 = m12;
		cur21 = m21;
		cur22 = m22;
		cur_prob = EXACT_TEST_BIAS;
		while (cur11 > 0.5) {
			cur12 += 1;
			cur21 += 1;
			cur_prob *= (cur11 * cur22) / (cur12 * cur21);
			cur11 -= 1;
			cur22 -= 1;
			preaddp = left_prob;
			left_prob += cur_prob;
			if (left_prob <= preaddp) {
				break;
			}
		}
		return left_prob / (left_prob + right_prob);
	}
}

