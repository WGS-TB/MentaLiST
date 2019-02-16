using Distributed
@everywhere using JuMP, Gurobi
@everywhere function kmer_coverage_ilp(locus, kmer_class, allele_ids, coverages)
  # m = Model(solver=GurobiSolver(Presolve=0))
  # one minute at most per locus, gap=0.1, no need to have optimal per se;
  # One thread per solve, since every ILP is already being run in parallel.
  m = Model(solver=GurobiSolver(OutputFlag=0, MIPGap=0.1, TimeLimit=1200, Threads=1))
  n_kmer = length(kmer_class)
  n_alleles = length(allele_ids)
  # kmer decision variable:
  @variable(m, x[1:n_kmer], Bin)
  @variable(m, c[1:n_alleles], Int)
  # coverage constraints:
  kmer_list = collect(keys(kmer_class))
  kmer_cover = DefaultDict{Int,Vector{Int}}(() -> Int[]) # allele -> list of kmers that cover it
  for (i, kmer) in enumerate(kmer_list)
    for al in kmer_class[kmer]
      push!(kmer_cover[al], i)
    end
  end
  # define coverage:
  for j in 1:n_alleles
    @constraint(m, c[j] == sum(x[i] for i in kmer_cover[j]))
  end
  # min coverage:
  for j in 1:n_alleles
    @constraint(m, c[j] >= coverages[j])
    # @constraint(m, c[j] == coverages[j])
  end
  # Forcing same coverage: (takes much longer, get many more kmers; better to store coverage)
  # Update: consider the cardinality of the kmers, to minimize it also:
  # @variable(m, max_card, Int) # max cardinality of the selected kmers
  # for (i, kmer) in enumerate(kmer_list)
  #   @constraint(m, max_card >= x[i] * length(kmer_class[kmer]))
  # end
  # minimize # of selected kmers:
  @objective(m, Min, sum(x[i] for i in 1:n_kmer))
  # @objective(m, Min, sum(x[i] for i in 1:n_kmer) + 100000*max_card)
  # TODO: e parameters works well, reduces variance but increases kmer number as expected.
  status = solve(m)
  if status == :Infeasible
    return nothing
  end
  # println("Objective value: ", getobjectivevalue(m))
  # for (i,kmer) in enumerate(kmer_list)
  #   if getvalue(x[i]) == 1
  #     println("$kmer : $(kmer_class[kmer])")
  #   end
  # end
  # println(join([Int(round(getvalue(c[j]),0)) for j in 1:n_alleles],", "))
  try
    kmers = [kmer_list[i] for i in 1:n_kmer if getvalue(x[i]) == 1]
    coverages = [Int16(round(getvalue(c[j]))) for j in 1:n_alleles]
    return kmers, coverages
  catch
    println("What happened??? -> LP status: $status LOCUS:$locus.")
    exit()
  end
  # return kmers, coverages, Int16(round(getvalue(max_card)))
end
