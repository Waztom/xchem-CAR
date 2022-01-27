import React, { useState, useEffect } from 'react';
import axios from 'axios';
import Spinner from 'react-bootstrap/Spinner';
import Image from 'react-bootstrap/Image';

const ProductImage = ({ reactionid }) => {
  // Use hooks instead of classes
  const [isLoading, setLoading] = useState(true);
  const [Product, setProduct] = useState([]);

  useEffect(() => {
    async function fetchData() {
      try {
        const request = await axios.get(`api/products?search=${reactionid}`);
        setProduct(request.data);
        setLoading(false);
      } catch (err) {
        console.log(err);
      }
    }
    fetchData();
  }, []);

  if (isLoading) {
    return (
      <Spinner animation="border" role="status">
        <span className="sr-only">Loading...</span>
      </Spinner>
    );
  }

  return Product.map((product) => <Image src={product.image} fluid />);
};

export default ProductImage;
