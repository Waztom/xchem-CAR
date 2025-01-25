export const formatPreferredVendorsString = vendors => {
  if (vendors.length === 1) {
    return vendors[0];
  }

  const lastVendor = vendors.at(-1);
  const remainingVendors = vendors.slice(0, -1);

  return `${remainingVendors.join(', ')} or ${lastVendor}`;
};
