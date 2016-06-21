#ifndef HEAPPQ_HH
#define HEAPPQ_HH

#include <vector>

#include "detail/PriorityQueue.hh"


template <typename ComparisonPolicy_fn>
class HeapPQ
{
private:

    typedef std::vector<Item> Tree;
    typedef std::size_t Value;

    ComparisonPolicy_fn cmp_fn;
    // Member Variables.
    std::vector<Index> vertexid_2_location_map;
    Tree tree;
    Index size;
    Index max_size;


    void exch(Index i, Index j)
    {
        Item t_pair;
        t_pair = tree[i];
        tree[i] = tree[j];
        tree[j] = t_pair;

        vertexid_2_location_map[tree[j].index] = j; //vertexid_2_location_map[j];
        vertexid_2_location_map[tree[i].index] = i;// t;
    }

    void fixUp(Tree &tree, Index k)
    {
        while ( (k > 1) && (cmp_fn)(tree[k/2] , tree[k]))
        {
            exch(k, k/2);
            k = k / 2;
        }
    }

    void fixDown(Tree &tree, Index k, const Index size)
    {
        while ( 2 * k <= size)
        {
            Index j = 2 * k;
            if ( j < size && cmp_fn(tree[j] , tree[j+1]))
            {
                j++;
            }
            if( !cmp_fn(tree[k], tree[j]))
            {
                break;
            }
            exch(k, j);
            k = j;
        }
    }

public:
    HeapPQ(Priority _max_priority, Index _max_size)
        :max_size(_max_size)
    {
        size = 0;
        tree.reserve(_max_size+1);
        vertexid_2_location_map.reserve(_max_size + 1);
    }


    bool empty() const
    {
        return size == 0;
    }

    void insert(Index vertex_id, Priority priority)
    {

        // std::cout << "p_queue.insert("<< vertex_id << "," << priority<< ")" << std::endl;
        tree[++size] = Item(vertex_id, priority);
        vertexid_2_location_map[vertex_id] = size;
        fixUp(tree,size);
    }

    Item top()
    {
        // std::cout << "Top (" << tree[1].first << "," << tree[1].second<< ")" << std::endl ;
        return tree[1];
    }

    void pop()
    {
        exch(1,size);
        fixDown(tree, 1, (size -1));
        size--;
    }

    Item get(Index vertex_id)
    {
        Index location = vertexid_2_location_map[vertex_id];
        return tree[location];
    }

    void increase(Index vertex_id)
    {
        Index location = vertexid_2_location_map[vertex_id];
        change(vertex_id,tree[location].priority  + 1);
    }

    void decrease(Index vertex_id)
    {
        Index location = vertexid_2_location_map[vertex_id];
        change(vertex_id,tree[location].priority  - 1);
    }

    void change(Index  vertex_id, Priority priority)
    {
        Index location = vertexid_2_location_map[vertex_id];
        // std::cout<< "Change (" << vertex_id << "," << priority <<
        // ")@" << location << std::endl;
        tree[location] = Item(vertex_id, priority);
        // std::cout << "------(" << vertex_id << "," << tree[location].second<< ")@ "<< location << std::endl;
        fixUp(tree, location);
        fixDown(tree, location,size);
        location = vertexid_2_location_map[vertex_id];
        // std::cout << "------(" << vertex_id << "," << tree[location].second<< ")@ "<< location << std::endl;
    }


    void remove(Index vertex_id)
    {
        Index location = vertexid_2_location_map[vertex_id];
        exch(location,size);
        size--;
        fixUp(tree,location);
        fixDown(tree,location,size);

    }

    void make_empty()
    {
        size = 0;
    }
};

#endif
