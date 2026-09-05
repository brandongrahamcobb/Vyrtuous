ALTER TABLE autoassign_roles ADD COLUMN id SERIAL;
ALTER TABLE autoassign_roles DROP CONSTRAINT autoassign_roles_pkey;
ALTER TABLE autoassign_roles ADD PRIMARY KEY (id);
