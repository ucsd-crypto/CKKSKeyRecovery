package main

import (
	"fmt"
	"math"
	"math/cmplx"
	"math/rand"
	"math/big"
	"time"
	"unsafe"
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/ckks"
)

func randomFloat(min, max float64) float64 {
	return min + rand.Float64()*(max-min)
}

func randomComplex(min, max float64) complex128 {
	return complex(randomFloat(min, max), randomFloat(min, max))
}

////////////////////////////////////////////////////////////////////////////////
// The linear key recovery attack, applied to homomorphically computed sigmoid
// of random complex numbers in [-8, 8], in the full packing mode
////////////////////////////////////////////////////////////////////////////////
func attack_sigmoid() {

	var err error

	// This demo is modified on top of the example code lattigo/sigmoid.
	// The following computation packs random 8192 float64 values in the range [-8, 8]
	// and approximates the function 1/(exp(-x) + 1) over the range [-8, 8].
	// Once the homomorphic computation is done, we decrypt and decode the results,
	// and then re-encode the resulting noisy complex numbers into a polynomial m',
	// and then try to compute the linear equation s' = c[a]^{-1} * (m' - c[b])

	rand.Seed(time.Now().UnixNano())

	// Scheme params
	params := ckks.DefaultParams[ckks.PN15QP827pq]

	encoder := ckks.NewEncoder(params)

	// Keys
	kgen := ckks.NewKeyGenerator(params)
	var sk *ckks.SecretKey
	var pk *ckks.PublicKey
	sk, pk = kgen.GenKeyPair()

	// Relinearization key
	var rlk *ckks.EvaluationKey
	rlk = kgen.GenRelinKey(sk)

	// Encryptor
	encryptor := ckks.NewEncryptorFromPk(params, pk)

	// Decryptor
	decryptor := ckks.NewDecryptor(params, sk)

	// Evaluator
	evaluator := ckks.NewEvaluator(params)

	// Values to encrypt
	values := make([]complex128, params.Slots())
	for i := range values {
		values[i] = complex(randomFloat(-8, 8), 0)
	}

	fmt.Printf("CKKS parameters: logN = %d, logQ = %d, levels = %d, scale= %f, sigma = %f \n",
		params.LogN(), params.LogQP(), params.MaxLevel()+1, params.Scale(), params.Sigma())
	fmt.Println()
	fmt.Printf("Values     : %6f %6f %6f %6f...\n",
		round(values[0]), round(values[1]), round(values[2]), round(values[3]))
	fmt.Println()

	// Plaintext creation and encoding process
	plaintext := encoder.EncodeNew(values, params.LogSlots())


	// Encryption process
	var ciphertext *ckks.Ciphertext
	ciphertext = encryptor.EncryptNew(plaintext)

	fmt.Println("Evaluation of the function 1/(exp(-x)+1) in the range [-8, 8] (degree of approximation: 32)")

	// Evaluation process
	// We approximate f(x) in the range [-8, 8] with a Chebyshev interpolant of 33 coefficients (degree 32).
	chebyapproximation := ckks.Approximate(f, -8, 8, 33)
	a := chebyapproximation.A()
	b := chebyapproximation.B()

	// Change of variable
	evaluator.MultByConst(ciphertext, 2/(b-a), ciphertext)
	evaluator.AddConst(ciphertext, (-a-b)/(b-a), ciphertext)
	evaluator.Rescale(ciphertext, params.Scale(), ciphertext)
	// We evaluate the interpolated Chebyshev interpolant on the ciphertext
	if ciphertext, err = evaluator.EvaluateCheby(ciphertext, chebyapproximation, rlk); err != nil {
		panic(err)
	}
	fmt.Println("Done... Consumed levels:", params.MaxLevel()-ciphertext.Level())

	// Computation of the reference values
	for i := range values {
		values[i] = f(values[i])
	}

	// Print results and comparison
	printDebug(params, ciphertext, values, decryptor, encoder)

	// Decrypt ciphertext
	level := ciphertext.Level()
	plaintextNoisy := decryptor.DecryptNew(ciphertext)

	valuesExact := make([]complex128, params.Slots())
	valuesExact = encoder.DecodeAndRound(plaintextNoisy, params.LogSlots(), plaintextNoisy.Scale()) // 1<<35

	valuesNoisy := make([]complex128, params.Slots())
	valuesNoisy = encoder.DecodeAndRound(plaintextNoisy, params.LogSlots(), 1<<35)
	if plaintextNoisy.Level() != level { panic("level doesn't match") }
	if !plaintextNoisy.IsNTT() { panic("decrypted plaintext should be in NTT form") }
	fmt.Printf("level = %d, scale = %f, slots = %d\n", level, math.Log2(plaintextNoisy.Scale()), params.Slots())
	fmt.Printf("added noise = %f\n", maxDiff(valuesExact, valuesNoisy))

	fmt.Println("Now try to recovery the secret key...\n")

	sc := 1<<5
	ciphertextScaled := evaluator.MultByConstNew(ciphertext, sc)
	fmt.Printf("Scale up...level = %d, scale = %f\n", ciphertextScaled.Level(), math.Log2(ciphertextScaled.Scale()))

	plaintextScaled := decryptor.DecryptNew(ciphertextScaled)
	valuesScaled := make([]complex128, params.Slots())
	valuesScaled = encoder.DecodeAndRound(plaintextScaled, params.LogSlots(), plaintextScaled.Scale()/(1<<5))
	fmt.Printf("valuesScaled = %f + %f I\n", real(valuesScaled[0]), imag(valuesScaled[1]))
	fmt.Printf("plaintextScaled level = %d, scale = %f\n", plaintextScaled.Level(), math.Log2(plaintextScaled.Scale()))

	{
		scInv := 1.0/float64(sc)
		for i := range valuesNoisy {
			// B := float64(1<<35)
			re := real(valuesScaled[i])*scInv
			im := imag(valuesScaled[i])*scInv
			valuesNoisy[i] = complex(re, im)
		}
		// evaluator.MultByConst(ciphertext, scInv, ciphertext)
		fmt.Println("Scale down... at levels:", ciphertext.Level(), ", scale = ", math.Log2(ciphertext.Scale()))
		// evaluator.Rescale(ciphertext, params.Scale(), ciphertext)
		// fmt.Println("Rescale... at levels:", ciphertext.Level(), ", scale = ", math.Log2(ciphertext.Scale()))
		// level = ciphertext.Level()
		// plaintextNoisy = decryptor.DecryptNew(ciphertext)
	}

	// Re-encode
	plaintextTest := ckks.NewPlaintext(params, level, ciphertext.Scale())
	encoder.Encode(plaintextTest, valuesNoisy, params.LogSlots()) // plaintextTest is not in NTT
	if plaintextTest.IsNTT() {
		panic("encoded plaintext should not be in NTT form")
	}

	// Create a ring data structure
	var ringQ *ring.Ring
	ringQ, _ = ring.NewRing(params.N(), params.Qi())

	// mEncode is the re-encoded polynomial, in NTT form
	mEncode := ringQ.NewPolyLvl(level)
	ringQ.NTTLvl(level, plaintextTest.Value()[0], mEncode) // mEncode in NTT



	mNoisy := ringQ.NewPolyLvl(level)						  // plaintextNoisy in non NTT form
	mError := ringQ.NewPolyLvl(level)						  // encoding error
	ringQ.InvNTTLvl(level, plaintextNoisy.Value()[0], mNoisy) // mNoisy is not in NTT
	ringQ.SubLvl(level, plaintextTest.Value()[0], mNoisy, mError)
	nm := inftyNorm(ringQ, mError, params.N())
	fmt.Printf("encoding error = %s\n", nm.String());
	printPoly(ringQ, mNoisy, level)
	printPoly(ringQ, plaintextTest.Value()[0], level)
	// Compute the linear equation m' - c[b] = c[a] * s
	if !ciphertext.IsNTT() { panic("ciphertext should be in NTT form") }
	rhs := ringQ.NewPolyLvl(level)
	ringQ.SubLvl(level, mEncode, ciphertext.Value()[0], rhs) // rhs is in NTT

	aInv := ringQ.NewPolyLvl(level)
	InvPolyNTT(ringQ, level, ciphertext.Value()[1], aInv) // aInv is in NTT

	sGuess := ringQ.NewPolyLvl(level)
	ringQ.MulCoeffsMontgomeryLvl(level, rhs, aInv, sGuess) // sGuess is in NTT
	// ringQ.InvNTTLvl(level, sGuess, sGuess)

	s := sk.Get().CopyNew()
	ringQ.InvMForm(s, s)
	// ringQ.InvNTT(s, s)
	if ringQ.EqualLvl(level, sGuess, s) {
		fmt.Printf("Found key!\n")
	} else {
		fmt.Printf("Failed\n")
	}
}

//////////////////////////////////////////////////
// Inversion in the cyclotomic ring
//////////////////////////////////////////////////

// modular exponential, taken from lattigo/utils/utils.go
func modExp(x, e, p uint64, bredParams []uint64) (result uint64) {
	result = 1
	for i := e; i > 0; i >>= 1 {
		if i&1 == 1 {
			result = ring.BRed(result, x, p, bredParams)
		}
		x = ring.BRed(x, x, p, bredParams)
	}
	return result
}

// compute p1^{-1} in the ring r, and return the result in p2
func InvPolyNTT(r *ring.Ring, level uint64, p1, p2 *ring.Poly) {
	for i := uint64(0); i < level+1; i++ {
		qi := r.Modulus[i]
		p1tmp, p2tmp := p1.Coeffs[i], p2.Coeffs[i]
		bredParams := r.BredParams[i]
		for j := uint64(0); j < r.N; j = j + 8 {
			x := (*[8]uint64)(unsafe.Pointer(&p1tmp[j]))
			y := (*[8]uint64)(unsafe.Pointer(&p2tmp[j]))

			y[0] = ring.MForm(modExp(x[0], qi-2, qi, bredParams), qi, bredParams)
			y[1] = ring.MForm(modExp(x[1], qi-2, qi, bredParams), qi, bredParams)
			y[2] = ring.MForm(modExp(x[2], qi-2, qi, bredParams), qi, bredParams)
			y[3] = ring.MForm(modExp(x[3], qi-2, qi, bredParams), qi, bredParams)
			y[4] = ring.MForm(modExp(x[4], qi-2, qi, bredParams), qi, bredParams)
			y[5] = ring.MForm(modExp(x[5], qi-2, qi, bredParams), qi, bredParams)
			y[6] = ring.MForm(modExp(x[6], qi-2, qi, bredParams), qi, bredParams)
			y[7] = ring.MForm(modExp(x[7], qi-2, qi, bredParams), qi, bredParams)
		}
	}
}


func inftyNorm(r *ring.Ring, p *ring.Poly, N uint64) *big.Int {
	max := new(big.Int)

	level := uint64(len(p.Coeffs)-1)

	bigintCoeffs := make([]*big.Int, N)
	r.PolyToBigint(p, bigintCoeffs)

	QBigInt := ring.NewUint(1)
	for i := range r.Modulus[:level+1]{
		QBigInt.Mul(QBigInt, ring.NewUint(r.Modulus[i]))
	}

	QHalfBigInt := new(big.Int)
	QHalfBigInt.Set(QBigInt)
	QHalfBigInt.Rsh(QBigInt, 1)

	// Centers and absolute values
	var sign int
	for i := range bigintCoeffs{
		sign = bigintCoeffs[i].Cmp(QHalfBigInt)
		if sign == 1 || sign == 0 {
			bigintCoeffs[i].Sub(bigintCoeffs[i], QBigInt)
			bigintCoeffs[i].Abs(bigintCoeffs[i])
		}
	}

	for i := uint64(0); i < r.N; i++ {
		if bigintCoeffs[i].Cmp(max) > 0 {
			max = bigintCoeffs[i]
		}
	}
	return max
}
	
func printPoly(r *ring.Ring, p *ring.Poly, level uint64) {
	for i := uint64(0); i < level+1; i++ {
		// qi := r.Modulus[i]
		fmt.Printf("[%d] = %d %d %d %d...\n", r.Modulus[i], p.Coeffs[i][0], p.Coeffs[i][1], p.Coeffs[i][2], p.Coeffs[i][3])
	}
}

func printComplex(values []complex128) {
	fmt.Printf("%6.10f %6.10f %6.10f %6.10f %6.10f %6.10f %6.10f %6.10f...\n",
		values[0], values[1], values[2], values[3], values[4], values[5], values[6], values[7])
}

func f(x complex128) complex128 {
	return 1 / (cmplx.Exp(-x) + 1)
}

func round(x complex128) complex128 {
	var factor float64
	factor = 100000000
	a := math.Round(real(x)*factor) / factor
	b := math.Round(imag(x)*factor) / factor
	return complex(a, b)
}

func printDebug(params *ckks.Parameters, ciphertext *ckks.Ciphertext, valuesWant []complex128, decryptor ckks.Decryptor, encoder ckks.Encoder) (valuesTest []complex128) {
	valuesTest = encoder.Decode(decryptor.DecryptNew(ciphertext), params.LogSlots())

	fmt.Println()
	fmt.Printf("Level: %d (logQ = %d)\n", ciphertext.Level(), params.LogQLvl(ciphertext.Level()))
	fmt.Printf("Scale: 2^%f\n", math.Log2(ciphertext.Scale()))
	fmt.Printf("ValuesTest: %6.10f %6.10f %6.10f %6.10f...\n", valuesTest[0], valuesTest[1], valuesTest[2], valuesTest[3])
	fmt.Printf("ValuesWant: %6.10f %6.10f %6.10f %6.10f...\n", valuesWant[0], valuesWant[1], valuesWant[2], valuesWant[3])
	fmt.Println()

	precStats := ckks.GetPrecisionStats(params, nil, nil, valuesWant, valuesTest, math.Exp2(53))
	fmt.Println(precStats.String())

	return
}

func main() {
	attack_sigmoid()
}

func maxDiff(a []complex128, b []complex128) (max float64) {
	max = 0
	for i := range a {
		d := a[i] - b[i]
		if math.Abs(real(d)) > max {
			max = math.Abs(real(d))
		}
		if math.Abs(imag(d)) > max {
			max = math.Abs(imag(d))
		}
	}
	return
}

