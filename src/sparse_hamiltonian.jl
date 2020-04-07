"""
Number of links for the boundary conditions.
"""
num_links(basis::AbstractSzbasis, boundary::BdryCond) = boundary == PBC ? basis.K : (basis.K - 1)

# NOTE: The number of links will depend on lattice dimension #

"""
Create a sparse Hamiltonian matrix for a PBC/OBC BH chain in 1D.

    H = -\\sum_{<i, j>} t_{i,j} (b_i^\\dagger b_j + b_i b_j^\\dagger) + (U/2) \\sum_i n_i (n_i - 1) - \\sum_i \\mu_i n_i
"""
function sparse_hamiltonian(basis::AbstractSzbasis, Ts::AbstractVector{Float64}, mus::AbstractVector{Float64}, U::Float64; boundary::BdryCond=PBC)
    
    end_site = num_links(basis, boundary)
    end_site = basis.K

    length(Ts) == end_site || error("Incorrect number of Ts: $(length(Ts)) != $(end_site)")
    length(mus) == basis.K || error("Incorrect number of mus: $(length(mus)) != $(basis.K)")

    # Surely, this is done this way because  
    # sparse matrices should be built at once
    rows = Int64[]
    cols = Int64[]
    elements = Float64[]
    
    # Define the linear size of the 2D lattice
    L = Int64(sqrt(basis.K))
    
    for (i, bra) in enumerate(basis)
        # Diagonal part
        Usum = 0
        musum = 0.0
        for j in 1:basis.K
            Usum += bra[j] * (bra[j]-1)
            musum += mus[j] * bra[j]
        end
        push!(rows, i)
        push!(cols, i)
        push!(elements, U*Usum/2.0 - musum)

        j_row = 0 # Keep track of site's row index
        j_col = 0 # Keep track of site's column index
        # Off-diagonal part
        for j in 1:end_site
            # Horizontal hopping 
            j_next = j+1
            j_col += 1 
            if j_col == L # Right edge PBC
                j_next = j - (L-1)
                #j_row += 1 # If we reached here, we're about to go to row below
                j_col = 0
            end
            # Vertical hopping
            j_down = j+L
            if j > L^2-L # Bottom edge PBC
                j_down = j - L*(L-1)
            end
            # Tunnel right, tunnel left,tunnel down, tunnel up
            for (site1, site2) in [(j, j_next), (j_next, j),(j, j_down), (j_down, j)]
                if bra[site1] > 0 # < ...,bra[site1]=n_site1,... | 
                    ket = copy(bra)
                    ket[site1] -= 1
                    ket[site2] += 1
                    if ket in basis # otherwise, zero element due to orthogonality
                        push!(rows, i)
                        push!(cols, serial_num(basis, ket))
                        push!(elements, -Ts[j] * sqrt(bra[site1]) * sqrt(bra[site2]+1))
                    end
                end
            end
        end
    end

    sparse(rows, cols, elements, length(basis), length(basis))
end

function sparse_hamiltonian(basis::AbstractSzbasis, Ts::AbstractVector{Float64}, U::Float64; boundary::BdryCond=PBC)
    sparse_hamiltonian(basis, Ts, zeros(basis.K), U, boundary=boundary)
end

function sparse_hamiltonian(basis::AbstractSzbasis, T::Float64, mus::AbstractVector{Float64}, U::Float64; boundary::BdryCond=PBC)
    sparse_hamiltonian(basis, fill(T, num_links(basis, boundary)), mus, U, boundary=boundary)
end

function sparse_hamiltonian(basis::AbstractSzbasis, T::Float64, U::Float64; boundary::BdryCond=PBC)
    sparse_hamiltonian(basis, fill(T, num_links(basis, boundary)), U, boundary=boundary)
end
