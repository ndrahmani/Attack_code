default("parisize", 64G);
default("timer", 1);
bit      = 1024 ; /* Primes bit size */
ns       = powers(2, 11);
verbose  = 0  ;  /* Set to 1 for printing messages, 0 to hide them */
NUM_TEST = 100  ;  /* Number of tests */
NUM_PRIM = 2  ;  /* Number of primes */

/* Function to safely read a binary file if it exists */
read_field(filename) = {
  my (result);
  iferr(
    result = read(filename),
    E,  /* If an error occurs (e.g., file not found) */
    result = 0;
  );
  result;
}

/* Function to write number field to a text file */
write_field(filename, field) = {
  write(filename, field);  /* Create and write the number field to the file */
}

{
  for (n_i = 1, #ns -1, 
    n = ns[n_i];
    filename = Str("cyclotomic_field_", n, ".txt");
    
    /* Attempt to read the number field from the file */
    K = read_field(filename);
    
    if (!K,
      
      /* File does not exist, compute the cyclotomic polynomial and number field */
      P = polcyclo(n);
      print1("Computing the ", n, "-th cyclotomic number field defined by the irreducible polynomial: ", P, " of degree: ", poldegree(P));
      K = nfinit(P);
      
      /* Store the number field in the file */
      write_field(filename, K);
      print("... Completed and stored.");
    , P = K.pol;
      print("Using Stored number field... \nThe ", n, "-th cyclotomic number field \ndefined by the irreducible polynomial: ", P, " of degree: ", poldegree(P));
    );

    if (verbose,
      print("/********************************");
      print("Generating starting primes since ");
      print("we don't have a quantum computer");
      print("*********************************/");
    );

    [S, F] = [0, 0];
    if (verbose,
      print("**************************");
      print("Number of tests: ", NUM_TEST);
      print("**************************");
    );

    t0 = getwalltime();
    t0_ = gettime();
    
    for (test = 1, NUM_TEST,
      /* Generate two random primes */
      kill(ps);
      ps = vector(NUM_PRIM, i, randomprime([2^(bit-1), 2^bit]));

      /* Decompose primes in the number field K */
      Ps = vector(#ps, i, idealprimedec(K, ps[i]));

      /* Generate a random ideal in K */
      A = idealmul(K, 1, 1);
      for (i = 1, #Ps,
        P_above_p_i = Ps[i][random([1, #Ps[i]])];
        P_exp_p_i =1 ;/* random([1, 32]);*/;
        A = idealmul(K, A,  P_above_p_i);/*idealpow(K, P_above_p_i, P_exp_p_i));*/
      );

      /* STARTING THE ATTACK */
      /************************/
      /* Calculate the norm of the ideal A and find prime divisors */
      t0 = getwalltime();
      normA = idealnorm(K, A);
      p_div_normA = [];
      for (i = 1, #ps,
        if (normA % ps[i] == 0,
          p_div_normA = concat(p_div_normA, ps[i])
        )
      );

      /* Verify all primes were used to construct A */
      if (verbose, print("Assuring that we got all the primes that we constructed A from: ", ps == p_div_normA));

      /* Decompose the primes that divide the norm of A */
      if (verbose, print("Assume now that we factored the norm of the ideal A then we decompose the primes and check for each prime divide A and get its valuation: "));
      ind_pow = [];
      for (i = 1, #p_div_normA,
        P_above_p = idealprimedec(K, p_div_normA[i]);
        for (j = 1, #P_above_p,
          e = idealval(K, A, P_above_p[j]);
          if (e != 0,
            ind_pow = concat(ind_pow, [[j, e]])
          )
        )
      );

      /* Reconstruct the ideal A */
      new_A = idealmul(K, 1, 1);
      for (i = 1, #Ps,
        [pi, e_pi] = ind_pow[i];
        P_above_p_i = Ps[i][pi];
        P_exp_p_i = e_pi;
        new_A = idealmul(K, new_A, idealpow(K, P_above_p_i, P_exp_p_i));
      );
      t1 = getwalltime();
      /* Check if reconstruction was successful */
      if (new_A == A,
        S += 1;
      ,
        F += 1;
      );
      kill(new_A);
      kill(Ps);
      kill(ps);
      kill(A)
    );

    t1_ = gettime();
    t1 = getwalltime();

    print("NUMBER OF TESTS  : ", NUM_TEST);
    print("SUCCESS RATE     : ", S/NUM_TEST * 100, " %");
    print("AVERAGE WALL TIME: ", strtime(ceil((t1-t0)/NUM_TEST)));
    print("AVERAGE CPU TIME : ", strtime(ceil((t1_-t0_)/NUM_TEST)));
    print("======================================================");
    print();
    print();
  );

  print("DONE");
}

\q
