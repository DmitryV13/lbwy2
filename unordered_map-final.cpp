

#include <iostream>
#include <valarray>
#include <cassert>
#include <algorithm>

#define XORRAND 3183507249
using std::function;
using std::cin;
using std::cout;
using std::string;

class invalid_object : public std::exception {
private:
    std::string msg;
public:
    invalid_object(const char* message) : msg(message) {}

    const char* what() const noexcept override {
        return msg.c_str();
    }
};

template<typename F, typename S>
class pair {
public:
    F first;
    S second;

    pair() = default;

    pair(F first_, S second_) : first(first_), second(second_) {}

    void put(F first_, S second_) {
        first = first_;
        second = second_;
    }

    friend bool operator!=(const pair &a, const pair &b) {
        return !(a == b);
    }

    friend bool operator==(const pair &a, const pair &b) {
        return (a.first == b.first && a.second == b.second);
    }
};

template<typename SN>
class dStack {
private:
    short int size;
    short int capacity;
    short int index;//for insert
    SN **doubleStackArray;
public:
    dStack() : index(0), size(0), capacity(5) {
        //doubleStackArray = (E *) malloc(5 * sizeof(E));
        doubleStackArray = new SN *[5];
    }

    dStack(const dStack &target) : size(target.size), capacity(target.capacity), index(target.index) {
        doubleStackArray = new SN *[capacity];
        for (int i = 0; i < target.size; ++i) {
            doubleStackArray[i] = target.doubleStackArray[i];
        }
    }

    ~dStack() {
        free(doubleStackArray);
        doubleStackArray = nullptr;
    }

    void alloc(int newCapacity) {
        SN **tmp = new SN *[newCapacity];
        for (int i = 0; i < capacity; i++) {
            tmp[i] = doubleStackArray[i];
        }
        delete[]doubleStackArray;
        doubleStackArray = tmp;
        capacity = static_cast<short>(newCapacity);
    }

    void push(SN *newObj) {
        //fool security
        if (newObj) {
            if (size >= capacity) alloc(capacity * 2);
            doubleStackArray[index++] = newObj;
            size++;
        }
    }

    SN *&peekFirst() {
        return doubleStackArray[index - 2];
    }

    SN *&peekLast() {
        return doubleStackArray[index - 1];
    }

    void popLast() {
        size = size - 1 > 0 ? size - 1 : 0;
        index = index - 1 > 0 ? index - 1 : 0;
    }

    bool empty() {
        //return index < 1;
        return size == 0;
    }

    bool hasOne() {
        return size == 1;
        //return (index == 1 && doubleStackArray);
    }
};

template<typename K, typename V>
class SpecialNode {
public:
    pair<const K, V> m_key_m_value;
    size_t m_hash;
    SpecialNode *next;
    SpecialNode *left;
    SpecialNode *right;

    explicit SpecialNode(const K k, V v, size_t h)
            : m_key_m_value(k,v), m_hash(h), next(nullptr), left(nullptr), right(nullptr) {
    }
};

template<typename K, typename V, typename SN=SpecialNode<K, V>, typename Stack=dStack<SN>>
class iterator {
private:
    SN *currentGeneral;
    SN *currentLocal;
    Stack dobStack;
public:
    explicit iterator(SN *&list) : currentGeneral(list) {
        if(currentGeneral!= nullptr) {
            findLastOnBranch(dobStack, currentGeneral);
            currentLocal = dobStack.peekLast();
        }
    }

    iterator(iterator &iter) = default;

    iterator() : currentLocal(nullptr), currentGeneral(nullptr), dobStack() {};

    ~iterator() = default;

    iterator &operator++() {
        if (dobStack.hasOne())
            dobStack.popLast();

        if (dobStack.empty()) {
            if(!currentGeneral){
                currentLocal = nullptr;
                return *this;
            }
            currentGeneral = currentGeneral->next;
            if (currentGeneral) {
                findLastOnBranch(dobStack, currentGeneral);
                currentLocal = dobStack.peekLast();
            } else {
                currentLocal = nullptr;
                return *this;
            }
        } else {
            if (dobStack.peekFirst()->right == dobStack.peekLast()) {
                currentLocal = dobStack.peekFirst();
                dobStack.popLast();
            } else {
                SN *tmp = dobStack.peekFirst();
                dobStack.popLast();
                if (tmp->right)
                    findLastOnBranch(dobStack, tmp->right);
                currentLocal = dobStack.peekLast();
            }
        }
        return *this;
    }

    iterator operator++(int) {
        iterator copy(*this);
        ++(*this);
        return copy;
    }

    // is searching from the parent node including it
    // if x is not NULL I call "findLastOnBranch" with x as second parameter
    // searching has snake form
    void findLastOnBranch(Stack &stack, SN *node) {
        if (!node->right && !node->left) {
            stack.push(node);
            return;
        }
        if (node->left) {
            stack.push(node);
            return findLastOnBranch(stack, node->left);
        }
        if (node->right) {
            stack.push(node);
            return findLastOnBranch(stack, node->right);
        }
    }

    pair<const K, V>& operator*() {
        return currentLocal->m_key_m_value;
    }

    pair<const K, V>* operator->() {
        return &currentLocal->m_key_m_value;
    }

    bool operator!=(const iterator &other) const {
        return !(*this == other);
    }

    bool operator==(const iterator &other) const {
        return currentLocal == other.currentLocal;
    }
};

template<typename K, typename V, typename SN=SpecialNode<K, V>, typename Stack=dStack<const SN>>
class citerator {
private:
    const SN *currentGeneral;
    const SN *currentLocal;
    Stack dobStack;
public:
    explicit citerator(const SN *&list) : currentGeneral(list) {
        if(currentGeneral!= nullptr) {
            findLastOnBranch(dobStack, currentGeneral);
            currentLocal = dobStack.peekLast();
        }
    }

    citerator(const citerator &iter)= default;

    citerator() : currentLocal(nullptr), currentGeneral(nullptr), dobStack() {};

    ~citerator() = default;

    citerator& operator++() {
        if (dobStack.hasOne())
            dobStack.popLast();

        if (dobStack.empty()) {
            if(!currentGeneral){
                currentLocal = nullptr;
                return *this;
            }
            currentGeneral = currentGeneral->next;//FIXME: initial map null
            if (currentGeneral) {
                findLastOnBranch(dobStack, currentGeneral);
                currentLocal = dobStack.peekLast();
            } else {
                currentLocal = nullptr;
                return *this;
            }
        } else {
            if (dobStack.peekFirst()->right == dobStack.peekLast()) {
                currentLocal = dobStack.peekFirst();
                dobStack.popLast();
            } else {
                const SN *tmp = dobStack.peekFirst();
                dobStack.popLast();
                if(tmp->right)
                    findLastOnBranch(dobStack, tmp->right);
                currentLocal = dobStack.peekLast();
            }
        }
        return *this;
    }

    citerator operator++(int) {
        citerator copy(*this);
        ++(*this);
        return copy;
    }

    void findLastOnBranch(Stack &stack, const SN *node) {
        if (!node->right && !node->left) {
            stack.push(node);
            return;
        }
        if (node->left) {
            stack.push(node);
            return findLastOnBranch(stack, node->left);
        }
        if (node->right) {
            stack.push(node);
            return findLastOnBranch(stack, node->right);
        }
    }

    const SN* getNode(){
        return currentLocal;
    }

    const pair<const K, V>& operator*() {
        return currentLocal->m_key_m_value;
    }

    const pair<const K, V>* operator->() {
        return &currentLocal->m_key_m_value;
    }

    bool operator!=(const citerator &other) const {
        return !(*this==other);
    }

    bool operator==(const citerator &other) const {
        return currentLocal == other.currentLocal;
    }
};

template<typename K, typename V, typename SN = SpecialNode<K, V>, typename map_iterator = iterator<K, V>, typename map_citerator = citerator<K, V>>
class abstract_data_t {
public:
    virtual ~abstract_data_t() = 0;

    [[nodiscard]] virtual size_t hashEvaluate(int key) const = 0;

    [[nodiscard]] virtual size_t hashEvaluate(string key) const = 0;

    virtual void put(const K key, V value) = 0;

    virtual void put(pair<const K, V> pair_ins) = 0;

    virtual bool empty() = 0;

    [[nodiscard]] virtual V at(const K &key) const = 0;

    virtual size_t length() = 0;

    virtual size_t capacity() = 0;

    virtual map_iterator &begin() = 0;

    virtual map_iterator &end() = 0;

    virtual map_citerator &cbegin() const = 0;

    virtual map_citerator &cend() const = 0;

    virtual map_iterator find(K key) = 0;

    virtual map_citerator cfind(K key) const = 0;

    virtual map_iterator insert(map_iterator iterator_insert, pair<const K, V> insert_pair) = 0;

    virtual map_iterator erase(map_iterator iterator_erase) = 0;

    virtual void putAtRearrange(SN *target, SN **nha, int mem, SN **newList) = 0;

    virtual void order(SN *node, int newMem, SN **nha, SN **nlon) = 0;

    virtual void add(SN *&root, K key, V value, size_t hash) = 0;

    virtual void add_hresize(SN *&root, K key, V value, size_t hash) = 0;

    virtual void copyNode(SN *initOb, SN *where) = 0;

    virtual void haresize(int mem) = 0;

    virtual void clearRoot(SN *root) = 0;

    virtual void clearList(SN *nodeList) = 0;

    virtual void clear() = 0;

    virtual bool contains(K key) = 0;

    virtual size_t count(K key) = 0;

    virtual void extend(const abstract_data_t<K,V>&) = 0;

    virtual size_t erase(const K& key) = 0;

    virtual V &operator[](K key) = 0;

    virtual bool operator==(abstract_data_t &b) = 0;

    virtual bool operator!=(abstract_data_t &b) = 0;
};

template<typename K, typename V, typename SN, typename map_iterator, typename map_citerator>
inline abstract_data_t<K, V, SN, map_iterator, map_citerator>::~abstract_data_t() = default;


template<typename K, typename V, typename SN = SpecialNode<K, V>, typename map_iterator = iterator<K, V>, typename map_citerator = citerator<K, V>>
class UnorderedMap : public abstract_data_t<K, V> {
private:
    int size;
    int m_capacity;
    const int load_factor = 1;

    SN **hashArray;
    SN *listOfNodes;
public:

    map_iterator &begin() override{
        map_iterator *tmp = new iterator<K, V>(listOfNodes);
        return *tmp;
    };

    map_iterator &end() override{
        return *new iterator<K, V>();
    };

    map_citerator &cbegin()const override{
        const SN* u=listOfNodes;
        map_citerator *tmp = new citerator<K, V>(u);
        return *tmp;
    };

    map_citerator &cend()const override{
        return *new citerator<K, V>();
    };

    map_citerator cfind(K key) const override{
        for (auto i = this->cbegin(); i != this->cend(); i++) {
            if (i->first == key)
                return i;
        }
        return this->cend();
    }

    map_iterator find(K key) override{
        for (auto i = this->begin(); i != this->end(); i++) {
            if (i->first == key)
                return i;
        }
        return this->end();
    }

    map_iterator insert(map_iterator iterator_insert, pair<const K, V> insert_pair) override{
        if(iterator_insert==this->end()){
            this->put(insert_pair);
            return this->find(insert_pair.first);
        }
        map_iterator tmp = find(insert_pair.first);
        if(tmp!=this->end()) {// there are values with such keys
            return this->end();
        }
        else{
            this->put(insert_pair);
            return this->find(insert_pair.first);
        }
    }

    map_iterator erase(map_iterator iterator_erase) override{
        for(auto i = this->begin(), tmp=this->begin(); i != end(); ++i, ++tmp){
            if(i==iterator_erase){
                ++tmp;
                this->erase(iterator_erase->first);
                return tmp;
            }
        }
        throw invalid_object("Transmitted iterator was invalid");
    }

    static size_t erase(abstract_data_t<K,V>& target, const K& key){
        target.erase(key);
    };

    static bool findByKey(map_iterator begin, const map_iterator &end, K key) {
        for (auto i = begin; i != end; i++) {
            if (i->first == key)
                return true;
        }
        return false;
    }

    UnorderedMap() : size(0), m_capacity(10), hashArray(nullptr), listOfNodes(nullptr) {
        haresize(m_capacity);
    };

    explicit UnorderedMap(int capacity) : size(0), m_capacity(capacity), hashArray(nullptr), listOfNodes(nullptr) {
        haresize(m_capacity);
    };

    UnorderedMap(const K arr[], int length) : size(0), m_capacity(length * 2), hashArray(nullptr),
                                              listOfNodes(nullptr) {
        haresize(m_capacity);//*2
        for (int i = 0; i < length; ++i) {
            put(i, arr[i]);
        }
    };

    UnorderedMap(std::initializer_list<pair<const K, V>> values): size(0), m_capacity(values.size()*2), hashArray(nullptr), listOfNodes(
            nullptr){
        haresize(m_capacity);
        for (auto i=values.begin(); i!=values.end(); ++i) {
            put(i->first, i->second);
        }
    }

    //UnorderedMap(std::initializer_list<std::pair<const K, V>> values): size(0), m_capacity(values.size()*2), hashArray(nullptr), listOfNodes(
    //        nullptr){
    //    haresize(m_capacity);
    //    if(values.begin())
    //        for (auto i=values.begin(); i!=values.end(); ++i) {
    //            put(i->first, i->second);
    //        }
    //}

    UnorderedMap(const UnorderedMap &initObject) : size(initObject.size), m_capacity(initObject.m_capacity),
                                                   hashArray(nullptr), listOfNodes(nullptr) {
        haresize(m_capacity);
        for (int i = 0; i < m_capacity; i++) hashArray[i] = nullptr;
        SN **initArray = initObject.hashArray;
        for (int i = m_capacity - 1; i >= 0; i--) {
            if (initArray[i]) {
                SN *tmp = new SpecialNode(initArray[i]->m_key_m_value.first, initArray[i]->m_key_m_value.second, initArray[i]->m_hash);
                hashArray[i] = tmp;
                copyNode(initArray[i], hashArray[i]);
            }
        }
        function<void(SN *, SN *&)> fillList = [&fillList, this](SN *fromList, SN *&listON) -> void {
            if (fromList->next) {
                fillList(fromList->next, listON);
            }
            hashArray[(fromList->m_hash) % m_capacity]->next = listON;
            listON = hashArray[(fromList->m_hash) % m_capacity];
        };
        fillList(initObject.listOfNodes, listOfNodes);
    };

    UnorderedMap(map_iterator begin, const map_iterator &end) : UnorderedMap() {
        for (auto i = begin; i != end; i++) {
            this->put(i->first, i->second);
        }
    }

    ~UnorderedMap() override {
        clear();
    };

    //for integers
    [[nodiscard]] size_t hashEvaluate(int key) const override {
        size_t result = floor(((int) sqrt(key) ^ XORRAND) * (10 - M_PI / (key / 2 + key)));
        if (result > 10000) {
            std::string numberStr = std::to_string(result);
            int del = (numberStr[5] + numberStr[4]) * 10 + (numberStr[4] + numberStr[1]);
            result = floor((double) result * del);
        }
        return result;
    };

    [[nodiscard]] size_t hashEvaluate(string key) const override {
        size_t result = 3 * key.size();
        for (int i = 0; i < key.size(); ++i) {
            result += key[i];
        }
        return result;
    }

    bool isCleared(){
        if(!hashArray && !listOfNodes && size==0 && m_capacity==0)
            return true;
        return false;
    }

    //insert value in map by the key
    void put(const K key, V value) override {
        if(isCleared()){
            m_capacity=10;
            haresize(m_capacity);
        }
        if (size / m_capacity >= load_factor) haresize(m_capacity * 2);
        size_t hash = hashEvaluate(key);
        int index = hash % m_capacity;
        if (hashArray[index] == nullptr) {
            SN *newEl = new SpecialNode(key, value, hash);
            size++;
            hashArray[index] = newEl;
            newEl->next = listOfNodes;
            listOfNodes = newEl;
        } else {
            SN *tmp = hashArray[index];
            add(tmp, key, value, hash);
        }
    };

    void put(pair<const K, V> pair_ins) override {
        if(isCleared()){
            m_capacity=10;
            haresize(m_capacity);
        }
        if (size / m_capacity >= load_factor) haresize(m_capacity * 2);
        size_t hash = hashEvaluate(pair_ins.first);
        int index = hash % m_capacity;
        if (hashArray[index] == nullptr) {
            SN *newEl = new SpecialNode(pair_ins.first, pair_ins.second, hash);
            size++;
            hashArray[index] = newEl;
            newEl->next = listOfNodes;
            listOfNodes = newEl;
        } else {
            SN *tmp = hashArray[index];
            add(tmp, pair_ins.first, pair_ins.second, hash);
        }
    };

    void extend(const abstract_data_t<K,V>& objectToCopy) override{
        if(isCleared()){
            m_capacity=10;
            haresize(m_capacity);
        }
        for(auto i =objectToCopy.cbegin(); i != objectToCopy.cend(); ++i){
            put(i->first, i->second);
        }
    }

    static bool is_equal(const UnorderedMap &object, UnorderedMap &target) {//meant const
        if (object.size == target.size && object.m_capacity == target.m_capacity) {
            SN *l1 = object.listOfNodes;
            while (l1) {
                if (!is_node_equal(l1, target)) return false;
                l1 = l1->next;
            }
            return true;
        }
        return false;
    };

    static int is_node_equal(SN *object, const UnorderedMap &target) {//meant const
        if (target.cfind(object->m_key_m_value.first) == target.cend()) return 0;
        int x = 1, y = 1;
        if (object->left) {
            x = is_node_equal(object->left, target);
        }
        if (object->right) {
            y = is_node_equal(object->right, target);
        }
        return x * y;
    }

    bool empty() override {
        return (!listOfNodes);
    };

    //swap containers' content
    void swap(UnorderedMap &other) {
        std::swap(m_capacity, other.m_capacity);
        std::swap(size, other.size);
        std::swap(hashArray, other.hashArray);
        std::swap(listOfNodes, other.listOfNodes);
    };

    //search by the key
    V at(const K &key) const override {
        int index = hashEvaluate(key) % m_capacity;
        function<SN *(SN *, const K)> find = [&find](SN *root, const K &key) -> SN * {
            if (root == nullptr) return nullptr;
            if (root->m_key_m_value.first == key) return root;
            if (key < root->m_key_m_value.first)
                return find(root->left, key);
            else if (key > root->m_key_m_value.first)
                return find(root->right, key);
            return nullptr;
        };
        SN *tmp = find(hashArray[index], key);
        if (!tmp) {
            throw std::out_of_range("Out of range");
        }
        return tmp->m_key_m_value.second;
    };//what does it return

    //how many inserted elements
    size_t length() override {
        return size;
    }

    //how many elements
    size_t capacity() override {
        return m_capacity;
    }

    //rearrangement of nodes in struct after hash array resize(changing hash array, list of nodes)
    void putAtRearrange(SN *target, SN **nha, int mem, SN **newList) override {
        size_t hash = target->m_hash;
        int index = hash % mem;
        if (nha[index] == nullptr) {
            nha[index] = target;
            target->next = *newList;
            *newList = target;
        } else {
            SN *tmp = nha[index];
            add_hresize(tmp, target->m_key_m_value.first, target->m_key_m_value.second, hash);//add
        }
    }

    //adding node in the right place in mini-tree after hash array resize
    void add(SN *&root, K key, V value, size_t hash) override {
        if (root == nullptr) {
            root = new SpecialNode(key, value, hash);
            size++;
        } else if (key < root->m_key_m_value.first) {
            add(root->left, key, value, hash);
        } else if (key > root->m_key_m_value.first) {
            add(root->right, key, value, hash);
        } else if (key == root->m_key_m_value.first) {
            root->m_key_m_value.second = value;
        }
    };

    void add_hresize(SN *&root, K key, V value, size_t hash) override{
        if (root == nullptr) {
            root = new SpecialNode(key, value, hash);
        } else if (key < root->m_key_m_value.first) {
            add_hresize(root->left, key, value, hash);
        } else if (key > root->m_key_m_value.first) {
            add_hresize(root->right, key, value, hash);
        } else if (key == root->m_key_m_value.first) {
            root->m_key_m_value.second = value;
        }
    };

    //copying mini-tree
    void copyNode(SN *initOb, SN *where) override {
        if (initOb->left != nullptr) {
            SN *objectForInsert = new SpecialNode(initOb->left->m_key_m_value.first, initOb->left->m_key_m_value.second,
                                                  initOb->left->m_hash);
            where->left = objectForInsert;
            copyNode(initOb->left, where->left);
        }
        if (initOb->right != nullptr) {
            SN *objectForInsert = new SpecialNode(initOb->right->m_key_m_value.first, initOb->right->m_key_m_value.second,
                                                  initOb->right->m_hash);
            where->right = objectForInsert;
            copyNode(initOb->right, where->right);
        }
    }

    //hash array resize
    void haresize(int mem) override {// mem should be > m_capacity
        SN **newHashArray = new SN *[mem];
        for (int i = 0; i < mem; i++) newHashArray[i] = nullptr;
        if (!hashArray) {
            hashArray = newHashArray;
            m_capacity = mem;
            return;
        }
        SN *newListOfNodes = nullptr;
        for (int i = 0; i < m_capacity; i++) {
            if (hashArray[i]) {
                order(hashArray[i], mem, newHashArray, &newListOfNodes);
            }
        }
        listOfNodes = newListOfNodes;
        m_capacity = mem;
        delete[]hashArray;
        hashArray = newHashArray;
    };

    size_t erase(const K& key) override{
        citerator target=this->cfind(key);
        if(target==this->cend())
            return 0;
        //after copying all children in array
        dStack<SN> stack;
        citerator helper(target);
        while(!helper.getNode()->next){
            ++helper;
        }
        hashArray[helper.getNode()->m_hash % m_capacity]= nullptr;
        saveRootElements(const_cast<SN*>(helper.getNode()), stack);
        eraseNodeFromListOfNodes(const_cast<SN*>(helper.getNode()));
        while(!stack.empty()){
            SN* tmp=stack.peekLast();
            stack.popLast();
            if(tmp==target.getNode()){
                hashArray[tmp->m_hash%m_capacity]= nullptr;
                delete tmp;
                tmp= nullptr;
            }else {
                putAtRearrange(tmp, hashArray, m_capacity, &listOfNodes);
            }
        }
        size--;
        return 1;
    }

    void saveRootElements(SN *node, dStack<SN>& stack){
        if (!node)
            return;
        if (node->left) {
            saveRootElements(node->left, stack);
            node->left = nullptr;
        }
        if (node->right) {
            saveRootElements(node->right, stack);
            node->right = nullptr;
        }
        stack.push(node);
    }

    void eraseNodeFromListOfNodes(SN* eraseTarget){
        SN* tmp=listOfNodes;
        if(tmp==eraseTarget){
            listOfNodes=listOfNodes->next;
            eraseTarget->next= nullptr;
            return;
        }
        while(tmp->next!=eraseTarget){
            tmp=tmp->next;
        }
        tmp->next=eraseTarget->next;
    }

    void order(SN *node, int newMem, SN **nha, SN **nlon) override {
        if (!node)
            return;
        if (node->left) {
            order(node->left, newMem, nha, nlon);
            node->left = nullptr;
        }
        if (node->right) {
            order(node->right, newMem, nha, nlon);
            node->right = nullptr;
        }
        putAtRearrange(node, nha, newMem, nlon);
    }

    //clearing mini-tree
    void clearRoot(SN *root) override {
        if (root == nullptr) return;
        if (root->left != nullptr) {
            clearRoot(root->left);
            root->left = nullptr;
        }
        if (root->right != nullptr) {
            clearRoot(root->right);
            root->right = nullptr;
        }
        delete root;
    }

    //clearing tree nodes
    void clearList(SN *nodeList) override {
        if (nodeList == nullptr) return;
        if (nodeList->next == nullptr) {
            clearRoot(nodeList);
        } else {
            clearList(nodeList->next);
            nodeList->next = nullptr;
            clearRoot(nodeList);
        }
    }

    void clear() override {
        clearList(listOfNodes);
        listOfNodes = nullptr;
        for (int i = 0; i < m_capacity; ++i) {
            if (hashArray[i]) {
                hashArray[i] = nullptr;
            }
        }
        delete[]hashArray;
        hashArray = nullptr;
        m_capacity = 0;
        size = 0;
    }

    bool contains(K key) override {
        at(key);
        return true;
    };

    size_t count(K key) override {
        at(key);
        return 1;
    };

    bool operator==(abstract_data_t<K, V> &b) override {//meant const
        return UnorderedMap::is_equal(*this, *dynamic_cast<UnorderedMap *>(&b));
    }

    bool operator!=(abstract_data_t<K, V> &b) override {
        return !(*this == b);
    }

    UnorderedMap &operator=(const UnorderedMap &initObj) {
        UnorderedMap tmp(initObj);
        swap(tmp);
        tmp.clear();
        return *this;
    };

    V &operator[](K key) override {
        if(isCleared()){
            m_capacity=10;
            haresize(m_capacity);
        }
        int index = hashEvaluate(key) % m_capacity;
        function<SN *(SN *, const K)> find = [&find](SN *root, const K &key) -> SN * {
            if (root == nullptr) return nullptr;
            if (root->m_key_m_value.first == key) return root;
            if (key < root->m_key_m_value.first)
                return find(root->left, key);
            else if (key > root->m_key_m_value.first)
                return find(root->right, key);
            return nullptr;
        };
        SN *tmp = find(hashArray[index], key);
        if (tmp == nullptr) {
            put(key, {});
            tmp = find(hashArray[index], key);
            return tmp->m_key_m_value.second;
        }
        return tmp->m_key_m_value.second;
    }

    friend std::ostream &operator<<(std::ostream &stream, const UnorderedMap<K,V> &a) {
        SN *tmp = a.listOfNodes;
        K *keys = static_cast<K *>(calloc(a.size, sizeof(K)));
        int i = 0;
        function<void(SN *)> findKeys = [&findKeys, &keys, &i](SN *root) -> void {
            if (root == nullptr) return;
            if (root->left != nullptr)
                findKeys(root->left);
            keys[i++] = root->m_key_m_value.first;
            if (root->right != nullptr)
                findKeys(root->right);
        };
        //finding keys
        for (; tmp; tmp = tmp->next) {
            findKeys(tmp);
        }
        //sorting keys
        for (i = 1; i < a.size; i++) {
            int targ = keys[i];
            int j = i - 1;
            while (j >= 0 && keys[j] > targ) {
                keys[j + 1] = keys[j];
                j = j - 1;
            }
            keys[j + 1] = targ;
        }
        //output
        for (i = 0; i < a.size; i++) {
            stream << a.at(keys[i]) << "\n";
        }
        free(keys);
        return stream;
    }

    friend std::istream &operator>>(std::istream &stream, UnorderedMap &a) {
        int tmp, i = 0;
        while (stream >> tmp) {
            a.put(i, tmp);
            i++;
        }
        return stream;
    }

};

template<typename K,typename V>
void print(const abstract_data_t<K,V>& a){
    std::cout<<*dynamic_cast<const UnorderedMap<K,V>*>(&a);
}

int main() {
    UnorderedMap<std::string, std::string> u ( {
            {"RED", "#FF0000"},
            {"GREEN", "#00FF00"},
            {"BLUE", "#0000FF"}});

    for (auto n : u)
        std::cout << n.first << ": " << n.second << " ";
    std::cout << std::endl;
    std::cout << std::endl;

    u["BLACK"] = "#000000";
    u["WHITE"] = "#FFFFFF";

    for (auto n : u)
        std::cout << n.first << ": " << n.second << " ";
}
