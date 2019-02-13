# Packages needed for parallel processing need to be loaded "twice":
using BioSequences
@everywhere using BioSequences

using DataStructures: DefaultDict
@everywhere using DataStructures: DefaultDict
using FastaIO
@everywhere using FastaIO
# using JuMP, Gurobi
# @everywhere using JuMP, Gurobi

@everywhere function twin(km)
    DNAKmer(reverse_complement(DNASequence(km)))
end

function contig_to_string(c)
  return "$(c[1])" * join(["$(x[end])" for x in c[2:end]])
end

@everywhere function get_contig(::Type{DNAKmer{k}}, kmers, km) where {k}
  contig_fw = get_contig_forward(DNAKmer{k}, kmers, km)
  contig_bw = get_contig_forward(DNAKmer{k}, kmers, twin(km))

  if km in neighbors(contig_fw[end])
    return contig_fw
  else
    return [[twin(x) for x in reverse(contig_bw)[1:end-1]] ; contig_fw]
  end
end

@everywhere function get_contig_forward(::Type{DNAKmer{k}}, kmers, km) where {k}
    c_fw = DNAKmer[km]
    while true
        # check if forward exists, and only 1:
        fw_neighbors = [kmer for kmer in neighbors(c_fw[end]) if canonical(kmer) in kmers]
        if length(fw_neighbors) != 1
          break
        end
        candidate = fw_neighbors[1]
        if candidate == km || candidate == twin(km) || candidate == twin(c_fw[end])
          break
        end
        bw_neighbors = [kmer for kmer in neighbors(twin(candidate)) if canonical(kmer) in kmers]
        if length(bw_neighbors) != 1
          break
        end
        push!(c_fw, candidate)
    end
    return c_fw
end

@everywhere function all_contigs(::Type{DNAKmer{k}}, kmers) where {k}
  done = Set{DNAKmer}()
  contigs = []
  for kmer in kmers
    if kmer in done
      continue
    end
    contig = get_contig(DNAKmer{k}, kmers, kmer)
      for c_kmer in contig
      push!(done, canonical(c_kmer))
    end
    push!(contigs, contig)
  end
  return contigs
end

### classifies each kmer with his colors (present alleles)
@everywhere function find_all_kmer_colors(::Type{DNAKmer{k}}, fastafile) where {k}
  # record = FASTASeqRecord{DNASequence}()
  record = FASTA.Record()
  kmer_class = DefaultDict{DNAKmer, Vector{Int16}}(() -> Int16[])
  allele_ids = Int16[]
  allele_idx::Int16 = 1
  n_kmers = Int[]
  seen = Set{DNAKmer}()
  # reader = open(FASTA.Reader{DNASequence}, fastafile)
  # reader = FASTA.Reader{DNASequence}, fastafile)
  reader = open(FASTA.Reader, fastafile)
  while !eof(reader)
    try
      read!(reader, record)
      seen = Set{DNAKmer}()
      for (pos, kmer) in each(DNAKmer{k}, FASTA.sequence(record))
        can_kmer = canonical(kmer)
        if !in(can_kmer, seen)
          push!(kmer_class[can_kmer], allele_idx)
          push!(seen, can_kmer)
        end
      end
    catch
      @error("CRITICAL ERROR: Error parsing file $fastafile, at record $(FASTA.identifier(record)), most likely some unkown characters. Please fix it and try again.")
      exit(-1)
    end
    push!(n_kmers, length(seen)) # number of unique kmers for this allele;
    # find the separator; will assume that if I see a "_", that's it, otherwise try "-";
    separator = in('_', FASTA.identifier(record)) ? "_" : "-"
    # update idx; the counter idx is incremental (1,2, ...) because we need the array sorted.
    # But this is not always sin the allele ordering, so we have to save the original id to restore it later;
    allele_id = parse(Int16,split(FASTA.identifier(record), separator)[end])
    push!(allele_ids, allele_id)
    allele_idx += 1
  end
  close(reader)
  return kmer_class, allele_ids, n_kmers
end

# @everywhere function kmer_coverage_ilp(locus, kmer_class, allele_ids, coverages)
#   # m = Model(solver=GurobiSolver(Presolve=0))
#   # one minute at most per locus, gap=0.1, no need to have optimal per se;
#   # One thread per solve, since every ILP is already being run in parallel.
#   m = Model(solver=GurobiSolver(OutputFlag=0, MIPGap=0.1, TimeLimit=1200, Threads=1))
#   n_kmer = length(kmer_class)
#   n_alleles = length(allele_ids)
#   # kmer decision variable:
#   @variable(m, x[1:n_kmer], Bin)
#   @variable(m, c[1:n_alleles], Int)
#   # coverage constraints:
#   kmer_list = collect(keys(kmer_class))
#   kmer_cover = DefaultDict{Int,Vector{Int}}(() -> Int[]) # allele -> list of kmers that cover it
#   for (i, kmer) in enumerate(kmer_list)
#     for al in kmer_class[kmer]
#       push!(kmer_cover[al], i)
#     end
#   end
#   # define coverage:
#   for j in 1:n_alleles
#     @constraint(m, c[j] == sum(x[i] for i in kmer_cover[j]))
#   end
#   # min coverage:
#   for j in 1:n_alleles
#     @constraint(m, c[j] >= coverages[j])
#     # @constraint(m, c[j] == coverages[j])
#   end
#   # Forcing same coverage: (takes much longer, get many more kmers; better to store coverage)
#   # Update: consider the cardinality of the kmers, to minimize it also:
#   # @variable(m, max_card, Int) # max cardinality of the selected kmers
#   # for (i, kmer) in enumerate(kmer_list)
#   #   @constraint(m, max_card >= x[i] * length(kmer_class[kmer]))
#   # end
#   # minimize # of selected kmers:
#   @objective(m, Min, sum(x[i] for i in 1:n_kmer))
#   # @objective(m, Min, sum(x[i] for i in 1:n_kmer) + 100000*max_card)
#   # TODO: e parameters works well, reduces variance but increases kmer number as expected.
#   status = solve(m)
#   if status == :Infeasible
#     return nothing
#   end
#   # println("Objective value: ", getobjectivevalue(m))
#   # for (i,kmer) in enumerate(kmer_list)
#   #   if getvalue(x[i]) == 1
#   #     println("$kmer : $(kmer_class[kmer])")
#   #   end
#   # end
#   # println(join([Int(round(getvalue(c[j]),0)) for j in 1:n_alleles],", "))
#   try
#     kmers = [kmer_list[i] for i in 1:n_kmer if getvalue(x[i]) == 1]
#     coverages = [Int16(round(getvalue(c[j]))) for j in 1:n_alleles]
#     return kmers, coverages
#   catch
#     println("What happened??? -> LP status: $status LOCUS:$locus.")
#     exit()
#   end
#   # return kmers, coverages, Int16(round(getvalue(max_card)))
# end


@everywhere function build_db_graph(::Type{DNAKmer{k}}, fastafile, coverage_p=1) where {k}
  # get kmers and colors (alleles)
  kmer_class, allele_ids, n_kmers = find_all_kmer_colors(DNAKmer{k}, fastafile)
  locus = splitext(basename(fastafile))[1]

  # Solve a coverage ILP to find the kmers:
  selected_kmers = DNAKmer[]
  actual_coverages = []
  # Coverage per allele:
  coverages = [ Int(round(n_k * coverage_p)) for n_k in n_kmers]
  # If coverage < 1, solve the ilp, otherwise just select all keys (all kmers)
  selected_kmers, actual_coverages = coverage_p < 1 ? kmer_coverage_ilp(locus, kmer_class, allele_ids, coverages) : (keys(kmer_class), coverages)
  # filter more kmers with de Bruijn graph contigs
  filtered_kmer_class, kmer_weights = filter_kmers_with_db_graph(DNAKmer{k}, kmer_class, selected_kmers)
  # debug log:
  sel_kmers, contig_kmers = length(selected_kmers), length(filtered_kmer_class)
  # println("$locus\t$coverage_p\tAll kmers\t$sel_kmers")
  # println("$locus\t$coverage_p\tdb Graph\t$contig_kmers")

  return (filtered_kmer_class, allele_ids, kmer_weights, actual_coverages)

end

@everywhere function filter_kmers_with_db_graph(::Type{DNAKmer{k}}, kmer_class, selected_kmers) where {k}
  # select 1 kmer per contig and set weight as length of contig; save kmers in filtered_kmer_class
  contig_list = all_contigs(DNAKmer{k}, selected_kmers)
  kmer_weights = Dict{UInt64, Int}() # Kmer is converted to UInt64, because DNAKmer apparently does not support write(), needed for the parallel pmap() call;
  filtered_kmer_class = Dict{UInt64, Vector{Int16}}()
  for contig in contig_list
      kmer = canonical(contig[1])
      kmer_uint = convert(UInt64,kmer)
      kmer_weights[kmer_uint] = length(contig)
      filtered_kmer_class[kmer_uint] = kmer_class[kmer]
  end
  return filtered_kmer_class, kmer_weights
end
